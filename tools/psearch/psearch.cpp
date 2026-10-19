// psearch.cpp - exact k-mer (k = 23) occurrence search over a refrel3 v1 cohort, three ways, with the same lookup:
//   fa     : every window of the source FASTA (.fa.gz, gunzip + FASTA parse), the reference truth for N = 50
//   v1     : every window of the fully decoded v1 archive (entropy decode + copy execution, streaming)
//   parse  : the reference scanned once; per assembly only the windows that are not inside one reference-contiguous run
//            (residual windows) are looked up; occurrences inside runs are transferred from the reference scan
//            (sorted runs x sorted reference occurrences)
// Window rule: canonical domain = upper-case bases of the contig records (lower case folded to upper); a window is
// searched only if all k bases are A/C/G/T (any other byte, N or IUPAC, excludes it); windows never cross contigs;
// strand + = window equals the pattern, - = window equals its reverse complement (k odd: no palindromes).
//   psearch gen <ref.fa> <archive dir> <names.txt> <seed> <out patterns.tsv>
//   psearch run <ref.fa> <archive dir> <names.txt> <patterns.tsv> <modes fa,v1,parse> <threads> [fa dir] [dump dir]
#include "pv_v1.h"
#include <zlib.h>
#include <atomic>
#include <chrono>
#include <map>
#include <unordered_set>
#include <array>
#include <mutex>
#include <random>
#include <sys/resource.h>
#include <thread>
using namespace pv;

static const int K = 23;
static const uint64_t KMASK = (1ull << (2 * K)) - 1;
static inline int b2(uint8_t c) { switch (c) { case 'A': return 0; case 'C': return 1; case 'G': return 2; case 'T': return 3; default: return -1; } }
static inline uint64_t rc_code(uint64_t x) { uint64_t r = 0; for (int i = 0; i < K; i++) { r = (r << 2) | (3 - (x & 3)); x >>= 2; } return r; }
static std::string decode_code(uint64_t x) { std::string s(K, 'A'); for (int i = K - 1; i >= 0; i--) { s[i] = "ACGT"[x & 3]; x >>= 2; } return s; }
static inline double now_s() { return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); }
static double peak_rss_gb() { struct rusage u; getrusage(RUSAGE_SELF, &u); return u.ru_maxrss / 1048576.0; }

// ---------------------------------------------------------------- lookup: open addressing + 2^20-bit prefilter
struct Entry { uint32_t pid; uint8_t strand; };
struct Table {
    std::vector<uint64_t> key; std::vector<uint32_t> first, count; std::vector<Entry> ent; std::vector<uint64_t> bits; uint64_t mask = 0;
    static inline uint64_t h1(uint64_t x) { return (x * 0x9E3779B97F4A7C15ull) >> 40; }
    static inline uint64_t h2(uint64_t x) { return ((x ^ (x >> 29)) * 0xBF58476D1CE4E5B9ull) >> 44; }
    void build(const std::vector<uint64_t>& fwd) {
        std::map<uint64_t, std::vector<Entry>> m;
        for (uint32_t p = 0; p < fwd.size(); p++) { m[fwd[p]].push_back({p, 0}); m[rc_code(fwd[p])].push_back({p, 1}); }
        size_t n = 1; while (n < 4 * m.size()) n <<= 1; mask = n - 1; key.assign(n, ~0ull); first.assign(n, 0); count.assign(n, 0); bits.assign((1u << 20) / 64, 0);
        for (auto& kv : m) { uint64_t i = h1(kv.first) & mask; while (key[i] != ~0ull) i = (i + 1) & mask;
            key[i] = kv.first; first[i] = (uint32_t)ent.size(); count[i] = (uint32_t)kv.second.size(); for (auto& e : kv.second) ent.push_back(e);
            const uint64_t b = h2(kv.first); bits[b >> 6] |= 1ull << (b & 63); }
    }
    inline int find(uint64_t x) const {   // slot or -1
        const uint64_t b = h2(x); if (!(bits[b >> 6] >> (b & 63) & 1)) return -1;
        uint64_t i = h1(x) & mask; for (;;) { if (key[i] == x) return (int)i; if (key[i] == ~0ull) return -1; i = (i + 1) & mask; } }
};

struct Hit { uint32_t contig; uint64_t pos; uint32_t pid; uint8_t strand;
    bool operator<(const Hit& o) const { return contig != o.contig ? contig < o.contig : pos != o.pos ? pos < o.pos : pid != o.pid ? pid < o.pid : strand < o.strand; }
    bool operator==(const Hit& o) const { return contig == o.contig && pos == o.pos && pid == o.pid && strand == o.strand; } };
static inline void emit(const Table& T, int slot, uint32_t contig, uint64_t pos, bool flip, std::vector<Hit>& out) {
    for (uint32_t j = 0; j < T.count[slot]; j++) { const Entry& e = T.ent[T.first[slot] + j]; out.push_back({contig, pos, e.pid, (uint8_t)(flip ? 1 - e.strand : e.strand)}); } }

// rolling scanner over a stream of bases with resets (contig starts, invalid bases)
struct Roller { uint64_t code = 0; int valid = 0; inline void reset() { valid = 0; code = 0; }
    inline bool push(uint8_t c) { const int v = b2(c); if (v < 0) { valid = 0; return false; } code = ((code << 2) | (uint64_t)v) & KMASK; if (valid < K) valid++; return valid == K; } };

static std::vector<std::string> read_names(const std::string& p) { std::vector<std::string> v; FILE* f = fopen(p.c_str(), "r"); if (!f) throw Err(p + ": cannot open");
    char b[4096]; while (fgets(b, sizeof b, f)) { std::string s(b); while (!s.empty() && (s.back() == '\n' || s.back() == '\r')) s.pop_back(); if (!s.empty()) v.push_back(s.substr(0, s.find('\t'))); } fclose(f); return v; }

struct Ref { std::vector<uint8_t> R; uint64_t n = 0; uint8_t sha[32]; };
static Ref load_ref(const char* p) { Ref r; r.R = slurp(p); r.n = fasta_bases_inplace(r.R.data(), r.R.size(), p); r.R.resize(r.n); r.R.shrink_to_fit(); rr_sha256(r.R.data(), r.n, r.sha); return r; }

// ---------------------------------------------------------------- per-assembly op walker (shared by gen and parse)
// Walks blocks in order: entropy decode (ops), copy execution into ext = [K-1 prefix | block], and reports, per op piece
// split at contig boundaries, the stream range, its kind and its run start (start of the merged reference run, or ~0).
struct Piece { uint64_t s, len; int kind; uint64_t run_start; uint64_t src; };
struct Stage { double tok = 0, exec = 0, scan = 0, transfer = 0, ref = 0, gz = 0, resolve = 0, hops = 0, rbytes = 0, maxdepth = 0; };

// ---------------------------------------------------------------- gen
static int cmd_gen(char** argv) {
    Ref ref = load_ref(argv[2]); const std::string dir = argv[3]; std::vector<std::string> names = read_names(argv[4]); std::mt19937_64 rng(strtoull(argv[5], 0, 10));
    const char* out = argv[6];
    struct Cand { std::string cls, seq, origin; };
    std::map<std::string, std::vector<Cand>> pool; std::map<std::string, uint64_t> seen;
    const std::map<std::string, size_t> want = {{"ref", 1000}, {"lit", 1000}, {"boundary", 200}, {"short_ins", 200}, {"self", 200}, {"contig_end", 200}, {"block_edge", 200}};
    auto offer = [&](const std::string& cls, const uint8_t* w, const std::string& origin) {
        for (int i = 0; i < K; i++) if (b2(w[i]) < 0) return;
        Cand c{cls, std::string((const char*)w, K), origin}; auto& v = pool[cls]; const uint64_t n = ++seen[cls]; const size_t cap = want.at(cls);
        if (v.size() < cap) v.push_back(c); else { const uint64_t j = rng() % n; if (j < cap) v[j] = c; } };
    // reference patterns
    while (pool["ref"].size() < 1000) { const uint64_t p = rng() % (ref.n - K); bool ok = true; for (int i = 0; i < K; i++) if (b2(ref.R[p + i]) < 0) { ok = false; break; }
        if (ok) pool["ref"].push_back({"ref", std::string((const char*)&ref.R[p], K), "ref:" + std::to_string(p)}); }
    seen["ref"] = 1000;
    std::vector<RrOp> ops(RR_MAXOPS); std::vector<uint8_t> lit(16384), ext(K - 1 + 16384);
    for (auto& nm : names) {
        const std::string p = dir + "/" + nm + ".q16k.rr3"; Archive X; open_v1(slurp(p), p, ref.sha, ref.n, X);
        std::vector<uint64_t> cstart, cend; for (auto& c : X.rec) { cstart.push_back(c.boff); cend.push_back(c.boff + c.len); }
        size_t ci = 0;
        // whole-assembly bytes for class sampling are not kept; windows are taken from ext (block + K-1 prefix)
        bool in = false; int rdir = 0; uint64_t rend = 0, nextref = 0, rstart = 0;
        uint64_t prev_kind_end = 0; int prev_kind = -1;
        for (uint64_t b = 0; b < X.nb; b++) {
            const uint64_t bs = b * X.Q; const uint32_t blen = (uint32_t)std::min<uint64_t>(X.Q, X.nbases - bs);
            if (b) memmove(ext.data(), ext.data() + X.Q, K - 1); else memset(ext.data(), 'N', K - 1);
            uint32_t o, l; R3Diag g; X.locate(b, &o, &l, &g);
            const int k = r3_decode_block(X.P + o, l, X.T, ref.R.data(), ref.n, blen, g, ops.data(), (uint32_t)ops.size(), lit.data(), (uint32_t)lit.size());
            if (k < 0) throw Err(nm + ": decode");
            exec_ops(ops.data(), k, lit.data(), ref.R.data(), ext.data() + K - 1);
            auto W = [&](uint64_t start) -> const uint8_t* {   // ext holds stream bytes [bs - (K-1), bs + blen)
                if ((int64_t)start < (int64_t)bs - (K - 1) || start + K > bs + blen) return nullptr;
                return ext.data() + (int64_t)start - ((int64_t)bs - (K - 1)); };
            auto contig_of = [&](uint64_t pos) { size_t c = 0; while (c < cend.size() && pos >= cend[c]) c++; return c; };
            auto within = [&](uint64_t start) { const size_t c = contig_of(start); return c < cend.size() && start + K <= cend[c]; };
            auto sample = [&](const std::string& cls, uint64_t start) { if (!within(start)) return; const uint8_t* w = W(start); if (w) offer(cls, w, nm + ":" + std::to_string(start)); };
            // block edge window
            if (b && bs >= K / 2) sample("block_edge", bs - K / 2);
            for (int j = 0; j < k; j++) { const RrOp& q = ops[j]; const uint64_t s = bs + q.dst, len = q.len;
                while (ci < cend.size() && s >= cend[ci]) ci++;
                const bool isref = q.kind == 1 || q.kind == 3;
                if (isref) { const int d = q.kind == 3; const uint64_t ps = d ? q.src : q.src;
                    const bool cont = in && rdir == d && rend == s && (d ? q.src + len == nextref : ps == nextref) && !(ci < cstart.size() && s == cstart[ci]);
                    if (!cont) { if (in || prev_kind >= 0) { if (s >= K / 2) sample("boundary", s - K / 2); } in = true; rdir = d; rstart = s; }
                    rend = s + len; nextref = d ? q.src : q.src + len; }
                else {
                    if (in && s >= K / 2) sample("boundary", s - K / 2);
                    in = false;
                    if (q.kind == 0) { if (len >= (uint64_t)K) sample("lit", s + (rng() % (len - K + 1))); else if (s + len >= (uint64_t)K) sample("short_ins", s + len - (K + len) / 2); }
                    if (q.kind == 2) { if (len >= (uint64_t)K) sample("self", s + (rng() % (len - K + 1))); else if (s >= (uint64_t)K) sample("self", s + len - K); }
                }
                prev_kind = q.kind; prev_kind_end = s + len; (void)prev_kind_end; (void)rstart; }
            for (size_t c = 0; c < cend.size(); c++) { if (cstart[c] >= bs && cstart[c] + K <= bs + blen) sample("contig_end", cstart[c]); if (cend[c] >= K && cend[c] - K >= bs && cend[c] <= bs + blen) sample("contig_end", cend[c] - K); }
            if (blen == X.Q) {} }
        fprintf(stderr, "gen %s done\n", nm.c_str());
    }
    FILE* f = fopen(out, "w"); fprintf(f, "pid\tclass\tseq\torigin\n"); uint32_t pid = 0; std::map<std::string, int> seenseq;
    for (auto& cls : {"ref", "lit", "boundary", "short_ins", "self", "contig_end", "block_edge"}) for (auto& c : pool[cls]) {
        if (seenseq.count(c.seq)) continue; seenseq[c.seq] = 1; fprintf(f, "%u\t%s\t%s\t%s\n", pid++, c.cls.c_str(), c.seq.c_str(), c.origin.c_str()); }
    fclose(f);
    for (auto& kv : pool) fprintf(stderr, "class %s: %zu (offered %llu)\n", kv.first.c_str(), kv.second.size(), (unsigned long long)seen[kv.first]);
    return 0;
}

// ---------------------------------------------------------------- searches
static void search_fa(const std::string& path, const Table& T, std::vector<Hit>& out, Stage& st) {
    const double t0 = now_s(); gzFile g = gzopen(path.c_str(), "rb"); if (!g) throw Err(path + ": cannot open"); gzbuffer(g, 1 << 20);
    std::vector<uint8_t> buf(1 << 22); Roller r; int64_t contig = -1; uint64_t pos = 0; bool header = false; double scan = 0;
    for (;;) { const int n = gzread(g, buf.data(), (unsigned)buf.size()); if (n <= 0) break; const double a = now_s();
        for (int i = 0; i < n; i++) { const uint8_t c = buf[i];
            if (header) { if (c == '\n') header = false; continue; }
            if (c == '>') { header = true; contig++; pos = 0; r.reset(); continue; }
            if (c == '\n' || c == '\r') continue;
            const uint8_t u = (c >= 'a' && c <= 'z') ? (uint8_t)(c - 32) : c;
            if (r.push(u)) { const int s = T.find(r.code); if (s >= 0) emit(T, s, (uint32_t)contig, pos + 1 - K, false, out); }
            pos++; }
        scan += now_s() - a; }
    gzclose(g); st.scan += scan; st.gz += now_s() - t0 - scan;
}

static void search_v1(const Archive& X, const Ref& ref, const Table& T, std::vector<Hit>& out, Stage& st) {
    std::vector<RrOp> ops(RR_MAXOPS); std::vector<uint8_t> lit(X.Q), blk(X.Q); Roller r; size_t ci = 0;
    std::vector<uint64_t> cend; for (auto& c : X.rec) cend.push_back(c.boff + c.len);
    while (ci < cend.size() && cend[ci] == 0) ci++;
    for (uint64_t b = 0; b < X.nb; b++) {
        const uint64_t bs = b * X.Q; const uint32_t blen = (uint32_t)std::min<uint64_t>(X.Q, X.nbases - bs);
        double a = now_s(); uint32_t o, l; R3Diag g; X.locate(b, &o, &l, &g);
        const int k = r3_decode_block(X.P + o, l, X.T, ref.R.data(), ref.n, blen, g, ops.data(), (uint32_t)ops.size(), lit.data(), (uint32_t)lit.size());
        if (k < 0) throw Err(X.path + ": decode"); double z = now_s(); st.tok += z - a; a = z;
        exec_ops(ops.data(), k, lit.data(), ref.R.data(), blk.data()); z = now_s(); st.exec += z - a; a = z;
        for (uint32_t i = 0; i < blen; i++) { const uint64_t p = bs + i;
            while (ci < cend.size() && p >= cend[ci]) { ci++; r.reset(); }
            if (r.push(blk[i])) { const int s = T.find(r.code); if (s >= 0) emit(T, s, (uint32_t)ci, p + 1 - K - X.rec[ci].boff, false, out); } }
        st.scan += now_s() - a; }
}

struct RefOcc { uint64_t pos; uint32_t slot; };
struct Run { uint64_t a, len, ref_lo; int dir; uint32_t contig; };

static void search_parse(const Archive& X, const Ref& ref, const Table& T, const std::vector<RefOcc>& occ, std::vector<Hit>& out, Stage& st) {
    std::vector<RrOp> ops(RR_MAXOPS); std::vector<uint8_t> lit(X.Q), ext(K - 1 + X.Q, 'N');
    std::vector<uint64_t> cstart, cend; for (auto& c : X.rec) { cstart.push_back(c.boff); cend.push_back(c.boff + c.len); }
    size_t ci = 0; while (ci < cend.size() && cend[ci] == 0) ci++;
    std::vector<Run> runs;
    bool in = false; int rdir = 0; uint64_t rstart = 0, rend = 0, nextref = 0, rlo = 0, rhi = 0;
    auto close_run = [&]() { if (in && rend - rstart >= (uint64_t)K) runs.push_back({rstart, rend - rstart, rdir ? rlo : rlo, rdir, (uint32_t)ci}); in = false; };
    // residual scan state
    uint64_t last_e = ~0ull; uint64_t code = 0; int valid = 0;
    for (uint64_t b = 0; b < X.nb; b++) {
        const uint64_t bs = b * X.Q; const uint32_t blen = (uint32_t)std::min<uint64_t>(X.Q, X.nbases - bs);
        double a = now_s(); uint32_t o, l; R3Diag g; X.locate(b, &o, &l, &g);
        const int k = r3_decode_block(X.P + o, l, X.T, ref.R.data(), ref.n, blen, g, ops.data(), (uint32_t)ops.size(), lit.data(), (uint32_t)lit.size());
        if (k < 0) throw Err(X.path + ": decode"); double z = now_s(); st.tok += z - a; a = z;
        if (b) memmove(ext.data(), ext.data() + X.Q, K - 1);
        exec_ops(ops.data(), k, lit.data(), ref.R.data(), ext.data() + K - 1); z = now_s(); st.exec += z - a; a = z;
        const uint8_t* base = ext.data() + K - 1;   // base[i] = stream byte bs + i, valid for i >= -(K-1)
        auto lookup_range = [&](uint64_t e0, uint64_t e1) {   // window ends e in [e0, e1), all inside contig ci
            for (uint64_t e = e0; e < e1; e++) {
                if (last_e == ~0ull || last_e + 1 != e) {   // restart from max(contig start, e - K + 1)
                    const uint64_t x0 = std::max<uint64_t>(cstart[ci], e + 1 >= (uint64_t)K ? e + 1 - K : 0);
                    code = 0; valid = 0;
                    for (uint64_t x = x0; x <= e; x++) { const int v = b2(base[(int64_t)x - (int64_t)bs]); if (v < 0) { valid = 0; code = 0; } else { code = ((code << 2) | (uint64_t)v) & KMASK; valid++; } }
                } else { const int v = b2(base[(int64_t)e - (int64_t)bs]); if (v < 0) { valid = 0; code = 0; } else { code = ((code << 2) | (uint64_t)v) & KMASK; if (valid < K) valid++; } }
                last_e = e;
                if (valid >= K) { const int s = T.find(code); if (s >= 0) emit(T, s, (uint32_t)ci, e + 1 - K - cstart[ci], false, out); } } };
        for (int j = 0; j < k; j++) { const RrOp& q = ops[j]; uint64_t s = bs + q.dst, len = q.len; uint64_t src = q.src;
            while (len) {
                while (ci < cend.size() && s >= cend[ci]) { close_run(); ci++; last_e = ~0ull; }
                const uint64_t piece = std::min<uint64_t>(len, cend[ci] - s);
                if (q.kind == 1 || q.kind == 3) {
                    const int d = q.kind == 3; const uint64_t ps = d ? src + len - piece : src;
                    const bool cont = in && rdir == d && rend == s && (d ? ps + piece == nextref : ps == nextref);
                    if (!cont) { close_run(); in = true; rdir = d; rstart = s; rlo = ps; rhi = ps + piece; }
                    rend = s + piece; nextref = d ? ps : ps + piece; if (d) rlo = ps; else rhi = ps + piece; (void)rhi;
                    // windows ending in this piece that start before the run start are residual
                    const uint64_t lim = std::min(s + piece, rstart + K - 1);
                    if (s < lim) lookup_range(s, lim);
                    if (!d) src += piece;
                } else { close_run(); lookup_range(s, s + piece); }
                s += piece; len -= piece; } }
        st.scan += now_s() - a; }
    while (ci < cend.size()) { close_run(); ci++; }
    // transfer: runs sorted by reference start x occurrences sorted by reference position
    const double a = now_s();
    std::sort(runs.begin(), runs.end(), [](const Run& x, const Run& y) { return x.ref_lo < y.ref_lo; });
    size_t lo = 0;
    static const bool brk = getenv("PSEARCH_BREAK_RC") != nullptr;   // negative control: drop reverse-complement transfers
    for (const Run& R_ : runs) {
        if (brk && R_.dir) continue;
        while (lo < occ.size() && occ[lo].pos < R_.ref_lo) lo++;
        const uint64_t last = R_.ref_lo + R_.len - K;
        for (size_t j = lo; j < occ.size() && occ[j].pos <= last; j++) {
            const uint64_t r = occ[j].pos;
            const uint64_t ap = R_.dir ? R_.a + (R_.ref_lo + R_.len - K - r) : R_.a + (r - R_.ref_lo);
            emit(T, occ[j].slot, R_.contig, ap - cstart[R_.contig], R_.dir != 0, out); } }
    st.transfer += now_s() - a;
}

// ---------------------------------------------------------------- origin: split hits by where they come from
// Same walk as search_parse (copies executed once per block), but every hit is attributed: "transfer" if it was carried
// over from the single reference scan through a reference-contiguous run, "residual" if it was found by looking at the
// assembly's own bytes (literals, self copies, windows straddling a run boundary). For residual hits the context
// +/- 20 bases (clamped to the contig) is hashed, so unique context sequences can be counted per pattern without
// keeping the strings. Not a timing path: the counters and the context hashing are extra work.
static const int PRE = 64;                                               // ext prefix: >= K-1+20 and > 20+K-1+20
struct Origin {
    std::vector<uint64_t> total, transfer, residual, asms;
    std::vector<std::vector<uint64_t>> ctx;                              // context hashes, deduplicated only at the end
    void init(size_t n) { total.assign(n, 0); transfer.assign(n, 0); residual.assign(n, 0); asms.assign(n, 0); ctx.resize(n); }
    void counters(const Origin& o) { for (size_t i = 0; i < total.size(); i++) { total[i] += o.total[i]; transfer[i] += o.transfer[i]; residual[i] += o.residual[i]; asms[i] += o.asms[i]; } }
};

static void search_origin(const Archive& X, const Ref& ref, const Table& T, const std::vector<RefOcc>& occ, std::vector<Hit>& out, Stage& st, Origin& O, std::vector<uint64_t>& res_local) {
    std::vector<RrOp> ops(RR_MAXOPS); std::vector<uint8_t> lit(X.Q), ext(PRE + X.Q, 'N');
    std::vector<uint64_t> cstart, cend; for (auto& c : X.rec) { cstart.push_back(c.boff); cend.push_back(c.boff + c.len); }
    size_t ci = 0; while (ci < cend.size() && cend[ci] == 0) ci++;
    std::vector<Run> runs; bool in = false; int rdir = 0; uint64_t rstart = 0, rend = 0, nextref = 0, rlo = 0;
    auto close_run = [&]() { if (in && rend - rstart >= (uint64_t)K) runs.push_back({rstart, rend - rstart, rlo, rdir, (uint32_t)ci}); in = false; };
    uint64_t last_e = ~0ull, code = 0; int valid = 0;
    struct Pend { uint64_t e; uint32_t ci, pid; }; std::vector<Pend> pend, pend_next;
    uint64_t bs = 0; uint32_t blen = 0;
    auto hash_ctx = [&](uint64_t e, uint32_t c) -> uint64_t {
        const uint64_t p = e + 1 - K;
        const uint64_t lo = std::max<uint64_t>(cstart[c], p >= 20 ? p - 20 : 0);
        const uint64_t hi = std::min<uint64_t>(std::min<uint64_t>(cend[c], e + 21), bs + blen);
        const uint8_t* base = ext.data() + PRE;                          // base[i] = stream byte bs + i
        return XXH3_64bits(base + ((int64_t)lo - (int64_t)bs), (size_t)(hi - lo)); };
    for (uint64_t b = 0; b < X.nb; b++) {
        bs = b * X.Q; blen = (uint32_t)std::min<uint64_t>(X.Q, X.nbases - bs);
        double a = now_s(); uint32_t o, l; R3Diag g; X.locate(b, &o, &l, &g);
        const int k = r3_decode_block(X.P + o, l, X.T, ref.R.data(), ref.n, blen, g, ops.data(), (uint32_t)ops.size(), lit.data(), (uint32_t)lit.size());
        if (k < 0) throw Err(X.path + ": decode"); double z = now_s(); st.tok += z - a; a = z;
        if (b) memmove(ext.data(), ext.data() + X.Q, PRE);
        exec_ops(ops.data(), k, lit.data(), ref.R.data(), ext.data() + PRE); z = now_s(); st.exec += z - a; a = z;
        for (const Pend& q : pend) O.ctx[q.pid].push_back(hash_ctx(q.e, q.ci));
        pend.clear();
        const uint8_t* base = ext.data() + PRE;
        auto emit_residual = [&](int slot, uint64_t e) {
            for (uint32_t j = 0; j < T.count[slot]; j++) { const Entry& en = T.ent[T.first[slot] + j];
                out.push_back({(uint32_t)ci, e + 1 - K - cstart[ci], en.pid, en.strand});
                O.total[en.pid]++; O.residual[en.pid]++; res_local[en.pid]++;
                if (e + 21 <= bs + blen) O.ctx[en.pid].push_back(hash_ctx(e, (uint32_t)ci)); else pend_next.push_back({e, (uint32_t)ci, en.pid}); } };
        auto lookup_range = [&](uint64_t e0, uint64_t e1) {
            for (uint64_t e = e0; e < e1; e++) {
                if (last_e == ~0ull || last_e + 1 != e) {
                    const uint64_t x0 = std::max<uint64_t>(cstart[ci], e + 1 >= (uint64_t)K ? e + 1 - K : 0); code = 0; valid = 0;
                    for (uint64_t x = x0; x <= e; x++) { const int v = b2(base[(int64_t)x - (int64_t)bs]); if (v < 0) { valid = 0; code = 0; } else { code = ((code << 2) | (uint64_t)v) & KMASK; valid++; } }
                } else { const int v = b2(base[(int64_t)e - (int64_t)bs]); if (v < 0) { valid = 0; code = 0; } else { code = ((code << 2) | (uint64_t)v) & KMASK; if (valid < K) valid++; } }
                last_e = e;
                if (valid >= K) { const int sl = T.find(code); if (sl >= 0) emit_residual(sl, e); } } };
        for (int j = 0; j < k; j++) { const RrOp& q = ops[j]; uint64_t s = bs + q.dst, len = q.len; uint64_t src = q.src;
            while (len) {
                while (ci < cend.size() && s >= cend[ci]) { close_run(); ci++; last_e = ~0ull; }
                const uint64_t piece = std::min<uint64_t>(len, cend[ci] - s);
                if (q.kind == 1 || q.kind == 3) {
                    const int d = q.kind == 3; const uint64_t ps = d ? src + len - piece : src;
                    const bool cont = in && rdir == d && rend == s && (d ? ps + piece == nextref : ps == nextref);
                    if (!cont) { close_run(); in = true; rdir = d; rstart = s; rlo = ps; }
                    rend = s + piece; nextref = d ? ps : ps + piece; if (d) rlo = ps;
                    const uint64_t lim = std::min(s + piece, rstart + K - 1); if (s < lim) lookup_range(s, lim);
                    if (!d) src += piece;
                } else { close_run(); lookup_range(s, s + piece); }
                s += piece; len -= piece; } }
        pend.swap(pend_next); pend_next.clear(); st.scan += now_s() - a; }
    for (const Pend& q : pend) O.ctx[q.pid].push_back(hash_ctx(q.e, q.ci));
    while (ci < cend.size()) { close_run(); ci++; }
    const double a = now_s();
    std::sort(runs.begin(), runs.end(), [](const Run& x, const Run& y) { return x.ref_lo < y.ref_lo; });
    size_t lo = 0;
    for (const Run& R_ : runs) { while (lo < occ.size() && occ[lo].pos < R_.ref_lo) lo++; const uint64_t last = R_.ref_lo + R_.len - K;
        for (size_t j = lo; j < occ.size() && occ[j].pos <= last; j++) { const uint64_t r = occ[j].pos;
            const uint64_t ap = R_.dir ? R_.a + (R_.ref_lo + R_.len - K - r) : R_.a + (r - R_.ref_lo);
            const int slot = (int)occ[j].slot;
            for (uint32_t t = 0; t < T.count[slot]; t++) { const Entry& en = T.ent[T.first[slot] + t];
                out.push_back({R_.contig, ap - cstart[R_.contig], en.pid, (uint8_t)(R_.dir ? 1 - en.strand : en.strand)});
                O.total[en.pid]++; O.transfer[en.pid]++; } } }
    st.transfer += now_s() - a;
}

// ---------------------------------------------------------------- parse2: residual windows without executing copies
// A byte of the stream is resolved on demand from the op that produced it: a literal from the block's literal bytes,
// a reference copy from the reference at the copy's coordinates (reverse complement through rr_comp), a self copy by
// following it to the earlier position of the same block (the source is always strictly earlier, so this terminates).
// Only the bytes of residual windows are resolved; the block is never materialised. The last K-1 bytes of each block
// are resolved into a tail buffer so that a window spanning the block boundary has its context.
struct OpIdx {
    mutable uint64_t hops = 0, rbytes = 0, maxdepth = 0;
    const RrOp* ops = nullptr; int n = 0; const uint8_t* lit = nullptr; const uint8_t* ref = nullptr;
    uint64_t bs = 0; const uint8_t* tail = nullptr;                      // tail[i] = stream byte bs - (K - 1) + i
    inline int find(uint32_t d) const { int lo = 0, hi = n - 1; while (lo < hi) { const int m = (lo + hi + 1) / 2; if (ops[m].dst <= d) lo = m; else hi = m - 1; } return lo; }
    // resolve with a monotone op cursor: no binary search while walking a contiguous range inside known ops
    inline uint8_t at_hint(int64_t x, int& j) const {
        if (x < (int64_t)bs) { const int64_t i = x - (int64_t)bs + (K - 1); return i >= 0 ? tail[i] : (uint8_t)'N'; }
        const uint32_t d = (uint32_t)(x - (int64_t)bs);
        while (j > 0 && ops[j].dst > d) j--;
        while (j + 1 < n && ops[j + 1].dst <= d) j++;
        const RrOp& q = ops[j]; const uint32_t off = d - q.dst; rbytes++;
        if (q.kind == 0) return lit[q.src + off];
        if (q.kind == 1) return ref[q.src + off];
        if (q.kind == 3) return rr_comp(ref[q.src + q.len - 1 - off]);
        // self copy: the source is earlier in the same block and usually inside a reference copy (forward or reverse
        // complement), so following the chain ends at a literal or at the reference. Negative control: shift the source.
        static const bool brk_self = getenv("PSEARCH_BREAK_SELF") != nullptr;
        return at((int64_t)bs + (int64_t)q.src + off + (brk_self ? 1 : 0));
    }
    inline uint8_t at(int64_t x) const {
        if (x < (int64_t)bs) { const int64_t i = x - (int64_t)bs + (K - 1); return i >= 0 ? tail[i] : (uint8_t)'N'; }
        uint32_t d = (uint32_t)(x - (int64_t)bs);
        for (int guard = 0; guard < 1 << 20; guard++) {
            hops++; if ((uint64_t)guard + 1 > maxdepth) maxdepth = guard + 1;
            const RrOp& q = ops[find(d)]; const uint32_t off = d - q.dst;
            if (q.kind == 0) return lit[q.src + off];
            if (q.kind == 1) return ref[q.src + off];
            if (q.kind == 3) return rr_comp(ref[q.src + q.len - 1 - off]);
            d = (uint32_t)q.src + off;                                   // self copy: strictly earlier in this block
        }
        return 'N';
    }
};

static void search_parse2(const Archive& X, const Ref& ref, const Table& T, const std::vector<RefOcc>& occ, std::vector<Hit>& out, Stage& st) {
    static const bool brk = getenv("PSEARCH_BREAK_CTX") != nullptr;      // negative control: tail context shifted by one
    std::vector<RrOp> ops(RR_MAXOPS); std::vector<uint8_t> lit(X.Q);
    std::vector<uint8_t> tail(K - 1, 'N'), tail_next(K - 1, 'N');
    std::vector<uint64_t> cstart, cend; for (auto& c : X.rec) { cstart.push_back(c.boff); cend.push_back(c.boff + c.len); }
    size_t ci = 0; while (ci < cend.size() && cend[ci] == 0) ci++;
    std::vector<Run> runs; bool in = false; int rdir = 0; uint64_t rstart = 0, rend = 0, nextref = 0, rlo = 0;
    auto close_run = [&]() { if (in && rend - rstart >= (uint64_t)K) runs.push_back({rstart, rend - rstart, rlo, rdir, (uint32_t)ci}); in = false; };
    uint64_t last_e = ~0ull, code = 0; int valid = 0; OpIdx ix; ix.lit = lit.data(); ix.ref = ref.R.data();
    for (uint64_t b = 0; b < X.nb; b++) {
        const uint64_t bs = b * X.Q; const uint32_t blen = (uint32_t)std::min<uint64_t>(X.Q, X.nbases - bs);
        double a = now_s(); uint32_t o, l; R3Diag g; X.locate(b, &o, &l, &g);
        const int k = r3_decode_block(X.P + o, l, X.T, ref.R.data(), ref.n, blen, g, ops.data(), (uint32_t)ops.size(), lit.data(), (uint32_t)lit.size());
        if (k < 0) throw Err(X.path + ": decode"); double z = now_s(); st.tok += z - a; a = z;
        ix.ops = ops.data(); ix.n = k; ix.bs = bs; ix.tail = tail.data();
        int hint = 0;
        auto lookup_range = [&](uint64_t e0, uint64_t e1) {
            for (uint64_t e = e0; e < e1; e++) {
                if (last_e == ~0ull || last_e + 1 != e) {
                    const uint64_t x0 = std::max<uint64_t>(cstart[ci], e + 1 >= (uint64_t)K ? e + 1 - K : 0); code = 0; valid = 0;
                    int h = hint;
                    for (uint64_t x = x0; x <= e; x++) { const int v = b2(ix.at_hint((int64_t)x, h)); if (v < 0) { valid = 0; code = 0; } else { code = ((code << 2) | (uint64_t)v) & KMASK; valid++; } }
                    hint = h;
                } else { const int v = b2(ix.at_hint((int64_t)e, hint)); if (v < 0) { valid = 0; code = 0; } else { code = ((code << 2) | (uint64_t)v) & KMASK; if (valid < K) valid++; } }
                last_e = e;
                if (valid >= K) { const int s = T.find(code); if (s >= 0) emit(T, s, (uint32_t)ci, e + 1 - K - cstart[ci], false, out); } } };
        for (int j = 0; j < k; j++) { const RrOp& q = ops[j]; uint64_t s = bs + q.dst, len = q.len; uint64_t src = q.src;
            while (len) {
                while (ci < cend.size() && s >= cend[ci]) { close_run(); ci++; last_e = ~0ull; }
                const uint64_t piece = std::min<uint64_t>(len, cend[ci] - s);
                if (q.kind == 1 || q.kind == 3) {
                    const int d = q.kind == 3; const uint64_t ps = d ? src + len - piece : src;
                    const bool cont = in && rdir == d && rend == s && (d ? ps + piece == nextref : ps == nextref);
                    if (!cont) { close_run(); in = true; rdir = d; rstart = s; rlo = ps; }
                    rend = s + piece; nextref = d ? ps : ps + piece; if (d) rlo = ps;
                    const uint64_t lim = std::min(s + piece, rstart + K - 1); if (s < lim) lookup_range(s, lim);
                    if (!d) src += piece;
                } else { close_run(); lookup_range(s, s + piece); }
                s += piece; len -= piece; } }
        z = now_s(); st.resolve += z - a; a = z;
        { int h = k - 1; for (int i = 0; i < K - 1; i++) { const int64_t x = (int64_t)(bs + blen) - (K - 1) + i; tail_next[i] = x >= 0 ? ix.at_hint(x, h) : (uint8_t)'N'; } }
        if (brk) { for (int i = K - 2; i > 0; i--) tail_next[i] = tail_next[i - 1]; }
        tail.swap(tail_next); st.scan += now_s() - a; }
    while (ci < cend.size()) { close_run(); ci++; }
    st.hops += (double)ix.hops; st.rbytes += (double)ix.rbytes; st.maxdepth = std::max(st.maxdepth, (double)ix.maxdepth);
    const double a = now_s();
    std::sort(runs.begin(), runs.end(), [](const Run& x, const Run& y) { return x.ref_lo < y.ref_lo; });
    size_t lo = 0;
    for (const Run& R_ : runs) { while (lo < occ.size() && occ[lo].pos < R_.ref_lo) lo++; const uint64_t last = R_.ref_lo + R_.len - K;
        for (size_t j = lo; j < occ.size() && occ[j].pos <= last; j++) { const uint64_t r = occ[j].pos;
            const uint64_t ap = R_.dir ? R_.a + (R_.ref_lo + R_.len - K - r) : R_.a + (r - R_.ref_lo);
            emit(T, occ[j].slot, R_.contig, ap - cstart[R_.contig], R_.dir != 0, out); } }
    st.transfer += now_s() - a;
}

static std::vector<RefOcc> scan_reference(const Ref& ref, const Table& T, int threads) {
    std::vector<std::vector<RefOcc>> part(threads); const uint64_t n = ref.n, chunk = (n + threads - 1) / threads; std::vector<std::thread> th;
    for (int t = 0; t < threads; t++) th.emplace_back([&, t] { const uint64_t a = t * chunk, b = std::min(n, a + chunk); Roller r;
        for (uint64_t p = a >= (uint64_t)(K - 1) ? a - (K - 1) : 0; p < std::min<uint64_t>(n, b + K - 1); p++) { if (r.push(ref.R[p])) { const uint64_t st = p + 1 - K; if (st >= a && st < b) { const int s = T.find(r.code); if (s >= 0) part[t].push_back({st, (uint32_t)s}); } } } });
    for (auto& x : th) x.join(); std::vector<RefOcc> occ; for (auto& v : part) occ.insert(occ.end(), v.begin(), v.end());
    std::sort(occ.begin(), occ.end(), [](const RefOcc& x, const RefOcc& y) { return x.pos < y.pos; }); return occ;
}

static std::string digest(const std::vector<Hit>& h) { RrSha256 s; for (auto& x : h) { uint8_t b[17]; memcpy(b, &x.contig, 4); memcpy(b + 4, &x.pos, 8); memcpy(b + 12, &x.pid, 4); b[16] = x.strand; s.update(b, 17); } uint8_t d[32]; s.final(d); return hex(d, 32); }

static int cmd_run(int argc, char** argv) {
    const double t_all = now_s();
    Ref ref = load_ref(argv[2]); const std::string dir = argv[3]; std::vector<std::string> names = read_names(argv[4]);
    std::vector<uint64_t> fwd; std::vector<std::string> cls;
    { FILE* f = fopen(argv[5], "r"); char b[4096]; fgets(b, sizeof b, f); while (fgets(b, sizeof b, f)) { char c[64], s[64]; unsigned pid; if (sscanf(b, "%u\t%63s\t%63s", &pid, c, s) == 3) { uint64_t x = 0; for (int i = 0; i < K; i++) x = (x << 2) | (uint64_t)b2((uint8_t)s[i]); fwd.push_back(x); cls.push_back(c); } } fclose(f); }
    Table T; T.build(fwd);
    const std::string modes = argv[6]; const int threads = atoi(argv[7]); const std::string fadir = argc > 8 ? argv[8] : "";
    printf("PATTERNS\t%zu\ttable slots %zu\tentries %zu\treference bases %llu\tload %.1f s\n", fwd.size(), T.key.size(), T.ent.size(), (unsigned long long)ref.n, now_s() - t_all); fflush(stdout);
    const size_t N = names.size();
    struct Sum { uint64_t n = 0; std::string d; std::map<std::string, uint64_t> pc; };
    std::map<std::string, std::vector<Sum>> res; std::vector<RefOcc> occ_keep;
    const size_t NPAT = fwd.size();
    std::vector<Origin> Oth(threads); for (auto& o : Oth) o.init(NPAT);
    std::vector<uint64_t> rl_scratch(NPAT, 0);
    auto one = [&](const std::string& mode, size_t i, const std::vector<RefOcc>& occ, std::vector<Hit>& H, Stage& S, int tid, std::vector<uint64_t>& res_local) {
        if (mode == "fa") search_fa(fadir + "/" + names[i].substr(3) + ".fa.gz", T, H, S);
        else { const std::string p = dir + "/" + names[i] + ".q16k.rr3"; const double a = now_s(); Archive X; open_v1(slurp(p), p, ref.sha, ref.n, X); S.gz += now_s() - a;
            if (mode == "v1") search_v1(X, ref, T, H, S); else if (mode == "origin") { std::fill(res_local.begin(), res_local.end(), 0); search_origin(X, ref, T, occ, H, S, Oth[tid], res_local); for (size_t q = 0; q < res_local.size(); q++) if (res_local[q]) Oth[tid].asms[q]++; }
            else if (mode == "parse2") search_parse2(X, ref, T, occ, H, S); else search_parse(X, ref, T, occ, H, S); }
        const double a = now_s(); std::sort(H.begin(), H.end()); if (mode == "parse" || mode == "parse2" || mode == "origin") S.transfer += now_s() - a; else S.scan += now_s() - a; };
    auto pass = [&](const std::string& mode) {
        std::vector<Sum> R(N); std::vector<Stage> S(threads); std::vector<std::vector<uint64_t>> rl(threads, std::vector<uint64_t>(NPAT, 0));
        std::atomic<size_t> next{0}; std::vector<RefOcc> occ; double tref = 0;
        const double w0 = now_s();
        if (mode == "parse" || mode == "parse2" || mode == "origin") { const double a = now_s(); occ = scan_reference(ref, T, threads); tref = now_s() - a; }
        std::vector<std::thread> th;
        for (int t = 0; t < threads; t++) th.emplace_back([&, t] { for (;;) { const size_t i = next.fetch_add(1); if (i >= N) break;
            std::vector<Hit> H; one(mode, i, occ, H, S[t], t, rl[t]); R[i].n = H.size(); R[i].d = digest(H); for (auto& h : H) R[i].pc[cls[h.pid]]++; } });
        for (auto& x : th) x.join(); const double wall = now_s() - w0;
        Stage tot; for (auto& s : S) { tot.tok += s.tok; tot.exec += s.exec; tot.scan += s.scan; tot.transfer += s.transfer; tot.gz += s.gz; tot.resolve += s.resolve; tot.hops += s.hops; tot.rbytes += s.rbytes; tot.maxdepth = std::max(tot.maxdepth, s.maxdepth); }
        uint64_t hits = 0; for (auto& r : R) hits += r.n;
        printf("PASS\t%s\twall %.1f s\tthreads %d\tcpu-s: input %.1f\ttoken decode %.1f\tcopy exec %.1f\tresolve+scan %.1f\tscan %.1f\ttransfer+sort %.1f\treference scan wall %.1f s (%zu occurrences, %.2f GB)\thits %llu\tpeak RSS %.1f GB\tresolved bytes %.4g\tchain hops %.4g\tmax chain depth %.0f\n",
               mode.c_str(), wall, threads, tot.gz, tot.tok, tot.exec, tot.resolve, tot.scan, tot.transfer, tref, occ.size(), occ.size() * sizeof(RefOcc) / 1e9, (unsigned long long)hits, peak_rss_gb(), tot.rbytes, tot.hops, tot.maxdepth);
        fflush(stdout); if (mode == "parse" || mode == "parse2" || mode == "origin") occ_keep = std::move(occ); res[mode] = std::move(R); };
    size_t p0 = 0; std::vector<std::string> ml; while (p0 <= modes.size()) { const size_t q = modes.find(',', p0); ml.push_back(modes.substr(p0, q == std::string::npos ? std::string::npos : q - p0)); if (q == std::string::npos) break; p0 = q + 1; }
    for (auto& m : ml) pass(m);
    const std::string m0 = ml[0]; int fails = 0;
    for (size_t mi = 1; mi < ml.size(); mi++) { const std::string& m = ml[mi]; int shown = 0; uint64_t diff_asm = 0;
        for (size_t i = 0; i < N; i++) { if (res[m0][i].n == res[m][i].n && res[m0][i].d == res[m][i].d) continue; diff_asm++;
            if (shown >= 20) continue;
            std::vector<Hit> A, B; Stage s; std::vector<RefOcc> none; one(m0, i, occ_keep, A, s, 0, rl_scratch); one(m, i, occ_keep, B, s, 0, rl_scratch);
            std::vector<Hit> onlyA, onlyB; std::set_difference(A.begin(), A.end(), B.begin(), B.end(), std::back_inserter(onlyA)); std::set_difference(B.begin(), B.end(), A.begin(), A.end(), std::back_inserter(onlyB));
            for (auto& h : onlyA) if (shown < 20) { printf("DIFF\t%s only\t%s\tcontig %u\tpos %llu\tpid %u (%s)\tstrand %c\n", m0.c_str(), names[i].c_str(), h.contig, (unsigned long long)h.pos, h.pid, cls[h.pid].c_str(), h.strand ? '-' : '+'); shown++; }
            for (auto& h : onlyB) if (shown < 20) { printf("DIFF\t%s only\t%s\tcontig %u\tpos %llu\tpid %u (%s)\tstrand %c\n", m.c_str(), names[i].c_str(), h.contig, (unsigned long long)h.pos, h.pid, cls[h.pid].c_str(), h.strand ? '-' : '+'); shown++; } }
        printf("COMPARE\t%s vs %s\tassemblies differing %llu of %zu\t%s\n", m0.c_str(), m.c_str(), (unsigned long long)diff_asm, N, diff_asm ? "FAIL" : "PASS"); if (diff_asm) fails++; }
    { RrSha256 s; std::map<std::string, uint64_t> pc; for (size_t i = 0; i < N; i++) { s.update((const uint8_t*)res[m0][i].d.data(), res[m0][i].d.size()); for (auto& kv : res[m0][i].pc) pc[kv.first] += kv.second; }
      uint8_t d[32]; s.final(d); printf("DIGEST\t%s\tall assemblies\t%s\n", m0.c_str(), hex(d, 32).c_str()); for (auto& kv : pc) printf("HITS_BY_CLASS\t%s\t%llu\n", kv.first.c_str(), (unsigned long long)kv.second); }
    if (res.count("origin")) { Origin O; O.init(NPAT); for (auto& o : Oth) O.counters(o);
      std::vector<uint64_t> uniq(NPAT, 0); std::vector<uint64_t> buf;
      for (size_t i = 0; i < NPAT; i++) { buf.clear();
        for (auto& o : Oth) { buf.insert(buf.end(), o.ctx[i].begin(), o.ctx[i].end()); std::vector<uint64_t>().swap(o.ctx[i]); }
        std::sort(buf.begin(), buf.end()); uniq[i] = (uint64_t)(std::unique(buf.begin(), buf.end()) - buf.begin()); }
      const char* op = getenv("PSEARCH_ORIGIN_OUT"); FILE* f = fopen(op ? op : "/dev/null", "w");
      fprintf(f, "pid\tclass\ttotal\ttransfer\tresidual\tresidual_frac\tassemblies_with_residual\tunique_contexts\n");
      std::map<std::string, std::array<uint64_t, 5>> by;   // total, transfer, residual, patterns, patterns_with_residual
      uint64_t pat_res = 0, ctx_tot = 0;
      for (size_t i = 0; i < NPAT; i++) { const uint64_t t = O.total[i], tr = O.transfer[i], rr = O.residual[i];
        fprintf(f, "%zu\t%s\t%llu\t%llu\t%llu\t%.6f\t%llu\t%zu\n", i, cls[i].c_str(), (unsigned long long)t, (unsigned long long)tr, (unsigned long long)rr,
                t ? (double)rr / (double)t : 0.0, (unsigned long long)O.asms[i], (size_t)uniq[i]);
        auto& v = by[cls[i]]; v[0] += t; v[1] += tr; v[2] += rr; v[3]++; if (rr) { v[4]++; pat_res++; } ctx_tot += uniq[i]; }
      fclose(f);
      printf("ORIGIN_TOTAL\tpatterns %zu\twith residual hits %llu\tunique contexts (all patterns) %llu\n", NPAT, (unsigned long long)pat_res, (unsigned long long)ctx_tot);
      for (auto& kv : by) printf("ORIGIN_BY_CLASS\t%s\tpatterns %llu\twith residual %llu\ttotal %llu\ttransfer %llu\tresidual %llu\tresidual_frac %.6f\n",
            kv.first.c_str(), (unsigned long long)kv.second[3], (unsigned long long)kv.second[4], (unsigned long long)kv.second[0], (unsigned long long)kv.second[1],
            (unsigned long long)kv.second[2], kv.second[0] ? (double)kv.second[2] / (double)kv.second[0] : 0.0); }
    { FILE* f = fopen(getenv("PSEARCH_PER_ASM") ? getenv("PSEARCH_PER_ASM") : "/dev/null", "w"); if (f) { for (size_t i = 0; i < N; i++) { fprintf(f, "%s", names[i].c_str()); for (auto& m : ml) fprintf(f, "\t%s %llu %s", m.c_str(), (unsigned long long)res[m][i].n, res[m][i].d.c_str()); fprintf(f, "\n"); } fclose(f); } }
    printf("TOTAL\twall %.1f s\tpeak RSS %.1f GB\t%s\n", now_s() - t_all, peak_rss_gb(), fails ? "FAIL" : "PASS");
    return fails ? 5 : 0;
}

int main(int argc, char** argv) {
    try { if (argc >= 7 && !strcmp(argv[1], "gen")) return cmd_gen(argv); if (argc >= 8 && !strcmp(argv[1], "run")) return cmd_run(argc, argv); }
    catch (const std::exception& e) { fprintf(stderr, "%s\n", e.what()); return 2; }
    fprintf(stderr, "usage: see header\n"); return 1;
}
