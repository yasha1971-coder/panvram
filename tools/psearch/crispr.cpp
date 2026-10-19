// crispr.cpp - SpCas9 off-target site search (20-nt spacer + NGG PAM, up to M mismatches in the spacer, PAM exact)
// over a refrel3 v1 cohort, three ways with the same matcher (fa: source FASTA; v1: full decode; parse: reference scanned
// once + residual windows + transfer through reference-contiguous runs; see psearch.cpp for the window rule).
// Site at window start p (23 bases, all A/C/G/T, inside one contig):
//   strand +: w[0..19] vs spacer with <= M mismatches, w[20] any base, w[21..22] = GG
//   strand -: w[0..1] = CC, w[2] any base, w[3..22] vs reverse complement of the spacer with <= M mismatches
//   crispr guides <ref.fa> <archive dir> <names.txt> <patterns.tsv> <seed> <n_ref> <n_asm> <out guides.tsv>
//   crispr run <ref.fa> <archive dir> <names.txt> <guides.tsv> <M> <modes> <threads> [fa dir]
#include "pv_v1.h"
#include <zlib.h>
#include <atomic>
#include <chrono>
#include <map>
#include <set>
#include <random>
#include <sys/resource.h>
#include <thread>
using namespace pv;

static const int K = 23;
static const uint64_t KMASK = (1ull << (2 * K)) - 1, M40 = (1ull << 40) - 1, EVEN40 = 0x5555555555ull;
static inline int b2(uint8_t c) { switch (c) { case 'A': return 0; case 'C': return 1; case 'G': return 2; case 'T': return 3; default: return -1; } }
static inline uint64_t rc_n(uint64_t x, int n) { uint64_t r = 0; for (int i = 0; i < n; i++) { r = (r << 2) | (3 - (x & 3)); x >>= 2; } return r; }
static inline double now_s() { return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); }
static double peak_rss_gb() { struct rusage u; getrusage(RUSAGE_SELF, &u); return u.ru_maxrss / 1048576.0; }

struct Hit { uint32_t contig; uint64_t pos; uint32_t gid; uint8_t strand, mm;
    bool operator<(const Hit& o) const { return contig != o.contig ? contig < o.contig : pos != o.pos ? pos < o.pos : gid != o.gid ? gid < o.gid : strand < o.strand; } };
struct Matcher { std::vector<uint64_t> fw, rv; int M = 0;
    // emits (gid, strand, mm) for one 23-mer code; flip swaps the strand (reverse-complement transfer)
    template <class F> inline void match(uint64_t code, F&& f) const {
        if ((code & 0xF) == 0xA) { const uint64_t s = code >> 6; for (uint32_t g = 0; g < fw.size(); g++) { const uint64_t x = s ^ fw[g]; const int mm = __builtin_popcountll((x | (x >> 1)) & EVEN40); if (mm <= M) f(g, 0, mm); } }
        if ((code >> 42) == 5) { const uint64_t s = code & M40; for (uint32_t g = 0; g < rv.size(); g++) { const uint64_t x = s ^ rv[g]; const int mm = __builtin_popcountll((x | (x >> 1)) & EVEN40); if (mm <= M) f(g, 1, mm); } } } };

struct Roller { uint64_t code = 0; int valid = 0; inline void reset() { valid = 0; code = 0; }
    inline bool push(uint8_t c) { const int v = b2(c); if (v < 0) { valid = 0; return false; } code = ((code << 2) | (uint64_t)v) & KMASK; if (valid < K) valid++; return valid == K; } };
static std::vector<std::string> read_names(const std::string& p) { std::vector<std::string> v; FILE* f = fopen(p.c_str(), "r"); if (!f) throw Err(p + ": cannot open");
    char b[4096]; while (fgets(b, sizeof b, f)) { std::string s(b); while (!s.empty() && (s.back() == '\n' || s.back() == '\r')) s.pop_back(); if (!s.empty()) v.push_back(s.substr(0, s.find('\t'))); } fclose(f); return v; }
struct Ref { std::vector<uint8_t> R; uint64_t n = 0; uint8_t sha[32]; };
static Ref load_ref(const char* p) { Ref r; r.R = slurp(p); r.n = fasta_bases_inplace(r.R.data(), r.R.size(), p); r.R.resize(r.n); r.R.shrink_to_fit(); rr_sha256(r.R.data(), r.n, r.sha); return r; }
struct Stage { double tok = 0, exec = 0, scan = 0, transfer = 0, gz = 0; };

static void search_fa(const std::string& path, const Matcher& G, std::vector<Hit>& out, Stage& st) {
    const double t0 = now_s(); gzFile g = gzopen(path.c_str(), "rb"); if (!g) throw Err(path + ": cannot open"); gzbuffer(g, 1 << 20);
    std::vector<uint8_t> buf(1 << 22); Roller r; int64_t contig = -1; uint64_t pos = 0; bool header = false; double scan = 0;
    for (;;) { const int n = gzread(g, buf.data(), (unsigned)buf.size()); if (n <= 0) break; const double a = now_s();
        for (int i = 0; i < n; i++) { const uint8_t c = buf[i];
            if (header) { if (c == '\n') header = false; continue; }
            if (c == '>') { header = true; contig++; pos = 0; r.reset(); continue; }
            if (c == '\n' || c == '\r') continue;
            const uint8_t u = (c >= 'a' && c <= 'z') ? (uint8_t)(c - 32) : c;
            if (r.push(u)) G.match(r.code, [&](uint32_t gid, int s, int mm) { out.push_back({(uint32_t)contig, pos + 1 - K, gid, (uint8_t)s, (uint8_t)mm}); });
            pos++; }
        scan += now_s() - a; }
    gzclose(g); st.scan += scan; st.gz += now_s() - t0 - scan;
}
static void search_v1(const Archive& X, const Ref& ref, const Matcher& G, std::vector<Hit>& out, Stage& st) {
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
            if (r.push(blk[i])) G.match(r.code, [&](uint32_t gid, int s, int mm) { out.push_back({(uint32_t)ci, p + 1 - K - X.rec[ci].boff, gid, (uint8_t)s, (uint8_t)mm}); }); }
        st.scan += now_s() - a; }
}
struct RefSite { uint64_t pos; uint32_t gid; uint8_t strand, mm; };
struct Run { uint64_t a, len, ref_lo; int dir; uint32_t contig; };
static void search_parse(const Archive& X, const Ref& ref, const Matcher& G, const std::vector<RefSite>& occ, std::vector<Hit>& out, Stage& st) {
    std::vector<RrOp> ops(RR_MAXOPS); std::vector<uint8_t> lit(X.Q), ext(K - 1 + X.Q, 'N');
    std::vector<uint64_t> cstart, cend; for (auto& c : X.rec) { cstart.push_back(c.boff); cend.push_back(c.boff + c.len); }
    size_t ci = 0; while (ci < cend.size() && cend[ci] == 0) ci++;
    std::vector<Run> runs; bool in = false; int rdir = 0; uint64_t rstart = 0, rend = 0, nextref = 0, rlo = 0;
    auto close_run = [&]() { if (in && rend - rstart >= (uint64_t)K) runs.push_back({rstart, rend - rstart, rlo, rdir, (uint32_t)ci}); in = false; };
    uint64_t last_e = ~0ull, code = 0; int valid = 0;
    for (uint64_t b = 0; b < X.nb; b++) {
        const uint64_t bs = b * X.Q; const uint32_t blen = (uint32_t)std::min<uint64_t>(X.Q, X.nbases - bs);
        double a = now_s(); uint32_t o, l; R3Diag g; X.locate(b, &o, &l, &g);
        const int k = r3_decode_block(X.P + o, l, X.T, ref.R.data(), ref.n, blen, g, ops.data(), (uint32_t)ops.size(), lit.data(), (uint32_t)lit.size());
        if (k < 0) throw Err(X.path + ": decode"); double z = now_s(); st.tok += z - a; a = z;
        if (b) memmove(ext.data(), ext.data() + X.Q, K - 1);
        exec_ops(ops.data(), k, lit.data(), ref.R.data(), ext.data() + K - 1); z = now_s(); st.exec += z - a; a = z;
        const uint8_t* base = ext.data() + K - 1;
        auto lookup_range = [&](uint64_t e0, uint64_t e1) {
            for (uint64_t e = e0; e < e1; e++) {
                if (last_e == ~0ull || last_e + 1 != e) { const uint64_t x0 = std::max<uint64_t>(cstart[ci], e + 1 >= (uint64_t)K ? e + 1 - K : 0); code = 0; valid = 0;
                    for (uint64_t x = x0; x <= e; x++) { const int v = b2(base[(int64_t)x - (int64_t)bs]); if (v < 0) { valid = 0; code = 0; } else { code = ((code << 2) | (uint64_t)v) & KMASK; valid++; } } }
                else { const int v = b2(base[(int64_t)e - (int64_t)bs]); if (v < 0) { valid = 0; code = 0; } else { code = ((code << 2) | (uint64_t)v) & KMASK; if (valid < K) valid++; } }
                last_e = e;
                if (valid >= K) G.match(code, [&](uint32_t gid, int s, int mm) { out.push_back({(uint32_t)ci, e + 1 - K - cstart[ci], gid, (uint8_t)s, (uint8_t)mm}); }); } };
        for (int j = 0; j < k; j++) { const RrOp& q = ops[j]; uint64_t s = bs + q.dst, len = q.len; uint64_t src = q.src;
            while (len) {
                while (ci < cend.size() && s >= cend[ci]) { close_run(); ci++; last_e = ~0ull; }
                const uint64_t piece = std::min<uint64_t>(len, cend[ci] - s);
                if (q.kind == 1 || q.kind == 3) { const int d = q.kind == 3; const uint64_t ps = d ? src + len - piece : src;
                    const bool cont = in && rdir == d && rend == s && (d ? ps + piece == nextref : ps == nextref);
                    if (!cont) { close_run(); in = true; rdir = d; rstart = s; rlo = ps; }
                    rend = s + piece; nextref = d ? ps : ps + piece; if (d) rlo = ps;
                    const uint64_t lim = std::min(s + piece, rstart + K - 1); if (s < lim) lookup_range(s, lim);
                    if (!d) src += piece; }
                else { close_run(); lookup_range(s, s + piece); }
                s += piece; len -= piece; } }
        st.scan += now_s() - a; }
    while (ci < cend.size()) { close_run(); ci++; }
    const double a = now_s();
    std::sort(runs.begin(), runs.end(), [](const Run& x, const Run& y) { return x.ref_lo < y.ref_lo; });
    size_t lo = 0;
    static const bool brk = getenv("CRISPR_BREAK_MM") != nullptr;   // negative control: transferred mismatch counts off by one
    for (const Run& R_ : runs) { while (lo < occ.size() && occ[lo].pos < R_.ref_lo) lo++; const uint64_t last = R_.ref_lo + R_.len - K;
        for (size_t j = lo; j < occ.size() && occ[j].pos <= last; j++) { const uint64_t r = occ[j].pos;
            const uint64_t ap = R_.dir ? R_.a + (R_.ref_lo + R_.len - K - r) : R_.a + (r - R_.ref_lo);
            out.push_back({R_.contig, ap - cstart[R_.contig], occ[j].gid, (uint8_t)(R_.dir ? 1 - occ[j].strand : occ[j].strand), (uint8_t)(brk && occ[j].mm == 3 ? 4 : occ[j].mm)}); } }
    st.transfer += now_s() - a;
}
static std::vector<RefSite> scan_reference(const Ref& ref, const Matcher& G, int threads) {
    std::vector<std::vector<RefSite>> part(threads); const uint64_t n = ref.n, chunk = (n + threads - 1) / threads; std::vector<std::thread> th;
    for (int t = 0; t < threads; t++) th.emplace_back([&, t] { const uint64_t a = t * chunk, b = std::min(n, a + chunk); Roller r;
        for (uint64_t p = a >= (uint64_t)(K - 1) ? a - (K - 1) : 0; p < std::min<uint64_t>(n, b + K - 1); p++) if (r.push(ref.R[p])) { const uint64_t st = p + 1 - K;
            if (st >= a && st < b) G.match(r.code, [&](uint32_t gid, int s, int mm) { part[t].push_back({st, gid, (uint8_t)s, (uint8_t)mm}); }); } });
    for (auto& x : th) x.join(); std::vector<RefSite> occ; for (auto& v : part) occ.insert(occ.end(), v.begin(), v.end());
    std::sort(occ.begin(), occ.end(), [](const RefSite& x, const RefSite& y) { return x.pos < y.pos; }); return occ;
}
static std::string digest(const std::vector<Hit>& h) { RrSha256 s; for (auto& x : h) { uint8_t b[18]; memcpy(b, &x.contig, 4); memcpy(b + 4, &x.pos, 8); memcpy(b + 12, &x.gid, 4); b[16] = x.strand; b[17] = x.mm; s.update(b, 18); } uint8_t d[32]; s.final(d); return hex(d, 32); }

static int cmd_guides(char** argv) {
    Ref ref = load_ref(argv[2]); std::mt19937_64 rng(strtoull(argv[6], 0, 10)); const int nref = atoi(argv[7]), nasm = atoi(argv[8]); FILE* f = fopen(argv[9], "w");
    fprintf(f, "gid\tsource\tspacer\torigin\n"); int g = 0;
    while (g < nref) { const uint64_t p = rng() % (ref.n - K); bool ok = true; for (int i = 0; i < K; i++) if (b2(ref.R[p + i]) < 0) { ok = false; break; }
        if (!ok || ref.R[p + 21] != 'G' || ref.R[p + 22] != 'G') continue; fprintf(f, "%d\tref\t%.20s\tref:%llu\n", g++, (const char*)&ref.R[p], (unsigned long long)p); }
    FILE* pf = fopen(argv[5], "r"); char b[4096]; fgets(b, sizeof b, pf); int na = 0;
    while (na < nasm && fgets(b, sizeof b, pf)) { char c[64], s[64], o[512]; unsigned pid; if (sscanf(b, "%u\t%63s\t%63s\t%511s", &pid, c, s, o) == 4 && !strcmp(c, "lit")) { fprintf(f, "%d\tassembly_literal\t%.20s\t%s\n", g++, s, o); na++; } }
    fclose(pf); fclose(f); fprintf(stderr, "guides: %d reference, %d from assembly literals\n", nref, na); return 0;
}

static int cmd_run(int argc, char** argv) {
    const double t_all = now_s(); Ref ref = load_ref(argv[2]); const std::string dir = argv[3]; std::vector<std::string> names = read_names(argv[4]);
    Matcher G; G.M = atoi(argv[6]); std::vector<std::string> src;
    { FILE* f = fopen(argv[5], "r"); char b[4096]; fgets(b, sizeof b, f); while (fgets(b, sizeof b, f)) { char so[64], s[64]; unsigned gid; if (sscanf(b, "%u\t%63s\t%63s", &gid, so, s) == 3) { uint64_t x = 0; for (int i = 0; i < 20; i++) x = (x << 2) | (uint64_t)b2((uint8_t)s[i]); G.fw.push_back(x); G.rv.push_back(rc_n(x, 20)); src.push_back(so); } } fclose(f); }
    const std::string modes = argv[7]; const int threads = atoi(argv[8]); const std::string fadir = argc > 9 ? argv[9] : "";
    printf("GUIDES\t%zu\tmax mismatches %d\tPAM NGG\tload %.1f s\n", G.fw.size(), G.M, now_s() - t_all); fflush(stdout);
    const size_t N = names.size(); struct Sum { uint64_t n = 0; std::string d; uint64_t bymm[8] = {0}; };
    std::map<std::string, std::vector<Sum>> res; std::vector<RefSite> occ_keep;
    auto one = [&](const std::string& mode, size_t i, const std::vector<RefSite>& occ, std::vector<Hit>& H, Stage& S) {
        if (mode == "fa") search_fa(fadir + "/" + names[i].substr(3) + ".fa.gz", G, H, S);
        else { const std::string p = dir + "/" + names[i] + ".q16k.rr3"; const double a = now_s(); Archive X; open_v1(slurp(p), p, ref.sha, ref.n, X); S.gz += now_s() - a;
            if (mode == "v1") search_v1(X, ref, G, H, S); else search_parse(X, ref, G, occ, H, S); }
        const double a = now_s(); std::sort(H.begin(), H.end()); if (mode == "parse") S.transfer += now_s() - a; else S.scan += now_s() - a; };
    auto pass = [&](const std::string& mode) {
        std::vector<Sum> R(N); std::vector<Stage> S(threads); std::atomic<size_t> next{0}; std::vector<RefSite> occ; double tref = 0; const double w0 = now_s();
        if (mode == "parse") { const double a = now_s(); occ = scan_reference(ref, G, threads); tref = now_s() - a; }
        std::vector<std::thread> th;
        for (int t = 0; t < threads; t++) th.emplace_back([&, t] { for (;;) { const size_t i = next.fetch_add(1); if (i >= N) break;
            std::vector<Hit> H; one(mode, i, occ, H, S[t]); R[i].n = H.size(); R[i].d = digest(H); for (auto& h : H) R[i].bymm[h.mm]++; } });
        for (auto& x : th) x.join(); const double wall = now_s() - w0;
        Stage tot; for (auto& s : S) { tot.tok += s.tok; tot.exec += s.exec; tot.scan += s.scan; tot.transfer += s.transfer; tot.gz += s.gz; }
        uint64_t hits = 0; for (auto& r : R) hits += r.n;
        printf("PASS\t%s\twall %.1f s\tthreads %d\tcpu-s: input %.1f\ttoken decode %.1f\tcopy exec %.1f\tscan %.1f\ttransfer+sort %.1f\treference scan wall %.1f s (%zu sites, %.2f GB)\tsites %llu\tpeak RSS %.1f GB\n",
               mode.c_str(), wall, threads, tot.gz, tot.tok, tot.exec, tot.scan, tot.transfer, tref, occ.size(), occ.size() * sizeof(RefSite) / 1e9, (unsigned long long)hits, peak_rss_gb());
        fflush(stdout); if (mode == "parse") occ_keep = std::move(occ); res[mode] = std::move(R); };
    size_t p0 = 0; std::vector<std::string> ml; while (p0 <= modes.size()) { const size_t q = modes.find(',', p0); ml.push_back(modes.substr(p0, q == std::string::npos ? std::string::npos : q - p0)); if (q == std::string::npos) break; p0 = q + 1; }
    for (auto& m : ml) pass(m);
    const std::string m0 = ml[0]; int fails = 0;
    for (size_t mi = 1; mi < ml.size(); mi++) { const std::string& m = ml[mi]; int shown = 0; uint64_t diff_asm = 0;
        for (size_t i = 0; i < N; i++) { if (res[m0][i].n == res[m][i].n && res[m0][i].d == res[m][i].d) continue; diff_asm++; if (shown >= 20) continue;
            std::vector<Hit> A, B; Stage s; one(m0, i, occ_keep, A, s); one(m, i, occ_keep, B, s);
            auto eq = [](const Hit& x, const Hit& y) { return !(x < y) && !(y < x) && x.mm == y.mm; };
            size_t ia = 0, ib = 0; while ((ia < A.size() || ib < B.size()) && shown < 20) {
                if (ib == B.size() || (ia < A.size() && A[ia] < B[ib])) { const Hit& h = A[ia++]; printf("DIFF\t%s only\t%s\tcontig %u\tpos %llu\tguide %u (%s)\tstrand %c\tmm %d\n", m0.c_str(), names[i].c_str(), h.contig, (unsigned long long)h.pos, h.gid, src[h.gid].c_str(), h.strand ? '-' : '+', h.mm); shown++; }
                else if (ia == A.size() || B[ib] < A[ia]) { const Hit& h = B[ib++]; printf("DIFF\t%s only\t%s\tcontig %u\tpos %llu\tguide %u (%s)\tstrand %c\tmm %d\n", m.c_str(), names[i].c_str(), h.contig, (unsigned long long)h.pos, h.gid, src[h.gid].c_str(), h.strand ? '-' : '+', h.mm); shown++; }
                else { if (!eq(A[ia], B[ib])) { printf("DIFF\tmismatch count\t%s\tpos %llu\n", names[i].c_str(), (unsigned long long)A[ia].pos); shown++; } ia++; ib++; } } }
        printf("COMPARE\t%s vs %s\tassemblies differing %llu of %zu\t%s\n", m0.c_str(), m.c_str(), (unsigned long long)diff_asm, N, diff_asm ? "FAIL" : "PASS"); if (diff_asm) fails++; }
    { RrSha256 s; uint64_t bymm[8] = {0}; for (size_t i = 0; i < N; i++) { s.update((const uint8_t*)res[m0][i].d.data(), res[m0][i].d.size()); for (int k = 0; k < 8; k++) bymm[k] += res[m0][i].bymm[k]; }
      uint8_t d[32]; s.final(d); printf("DIGEST\t%s\tall assemblies\t%s\n", m0.c_str(), hex(d, 32).c_str()); for (int k = 0; k <= G.M; k++) printf("SITES_BY_MISMATCHES\t%d\t%llu\n", k, (unsigned long long)bymm[k]); }
    { FILE* f = fopen(getenv("CRISPR_PER_ASM") ? getenv("CRISPR_PER_ASM") : "/dev/null", "w"); if (f) { for (size_t i = 0; i < N; i++) { fprintf(f, "%s", names[i].c_str()); for (auto& m : ml) fprintf(f, "\t%s %llu %s", m.c_str(), (unsigned long long)res[m][i].n, res[m][i].d.c_str()); fprintf(f, "\n"); } fclose(f); } }
    printf("TOTAL\twall %.1f s\tpeak RSS %.1f GB\t%s\n", now_s() - t_all, peak_rss_gb(), fails ? "FAIL" : "PASS");
    return fails ? 5 : 0;
}
// residual: unique 23-nt sequences of the sites that are NOT transferred from the reference scan, per guide, in
// canonical form (lesser of the window and its reverse complement), with how many assemblies carry each and one
// example location. The transfer step is skipped: only the assembly's own residual windows are looked at.
static inline uint64_t rc23(uint64_t x) { return rc_n(x, K); }
static std::string code_str(uint64_t x) { std::string s(K, 'A'); for (int i = K - 1; i >= 0; i--) { s[i] = "ACGT"[x & 3]; x >>= 2; } return s; }
struct ResRec { uint64_t occ = 0, asms = 0; uint8_t mm = 0; uint32_t a = 0, c = 0; uint64_t pos = 0; uint8_t strand = 0; };
using ResMap = std::map<std::pair<uint32_t, uint64_t>, ResRec>;   // (guide, canonical code) -> record

static void search_residual(const Archive& X, const Ref& ref, const Matcher& G, uint32_t ai, ResMap& R, std::set<std::pair<uint32_t, uint64_t>>& seen_here) {
    std::vector<RrOp> ops(RR_MAXOPS); std::vector<uint8_t> lit(X.Q), ext(K - 1 + X.Q, 'N');
    std::vector<uint64_t> cstart, cend; for (auto& c : X.rec) { cstart.push_back(c.boff); cend.push_back(c.boff + c.len); }
    size_t ci = 0; while (ci < cend.size() && cend[ci] == 0) ci++;
    bool in = false; int rdir = 0; uint64_t rstart = 0, rend = 0, nextref = 0;
    uint64_t last_e = ~0ull, code = 0; int valid = 0;
    for (uint64_t b = 0; b < X.nb; b++) {
        const uint64_t bs = b * X.Q; const uint32_t blen = (uint32_t)std::min<uint64_t>(X.Q, X.nbases - bs);
        uint32_t o, l; R3Diag g; X.locate(b, &o, &l, &g);
        const int k = r3_decode_block(X.P + o, l, X.T, ref.R.data(), ref.n, blen, g, ops.data(), (uint32_t)ops.size(), lit.data(), (uint32_t)lit.size());
        if (k < 0) throw Err(X.path + ": decode");
        if (b) memmove(ext.data(), ext.data() + X.Q, K - 1);
        exec_ops(ops.data(), k, lit.data(), ref.R.data(), ext.data() + K - 1);
        const uint8_t* base = ext.data() + K - 1;
        auto lookup_range = [&](uint64_t e0, uint64_t e1) {
            for (uint64_t e = e0; e < e1; e++) {
                if (last_e == ~0ull || last_e + 1 != e) {
                    const uint64_t x0 = std::max<uint64_t>(cstart[ci], e + 1 >= (uint64_t)K ? e + 1 - K : 0); code = 0; valid = 0;
                    for (uint64_t x = x0; x <= e; x++) { const int v = b2(base[(int64_t)x - (int64_t)bs]); if (v < 0) { valid = 0; code = 0; } else { code = ((code << 2) | (uint64_t)v) & KMASK; valid++; } }
                } else { const int v = b2(base[(int64_t)e - (int64_t)bs]); if (v < 0) { valid = 0; code = 0; } else { code = ((code << 2) | (uint64_t)v) & KMASK; if (valid < K) valid++; } }
                last_e = e;
                if (valid >= K) { const uint64_t cc = code; G.match(cc, [&](uint32_t gid, int strand, int mm) {
                    const uint64_t can = std::min(cc, rc23(cc)); const auto key = std::make_pair(gid, can);
                    ResRec& r = R[key]; if (!r.occ) { r.mm = (uint8_t)mm; r.a = ai; r.c = (uint32_t)ci; r.pos = e + 1 - K - cstart[ci]; r.strand = (uint8_t)strand; }
                    r.occ++; if (seen_here.insert(key).second) r.asms++; }); } } };
        for (int j = 0; j < k; j++) { const RrOp& q = ops[j]; uint64_t s = bs + q.dst, len = q.len; uint64_t src = q.src;
            while (len) {
                while (ci < cend.size() && s >= cend[ci]) { in = false; ci++; last_e = ~0ull; }
                const uint64_t piece = std::min<uint64_t>(len, cend[ci] - s);
                if (q.kind == 1 || q.kind == 3) { const int d = q.kind == 3; const uint64_t ps = d ? src + len - piece : src;
                    const bool cont = in && rdir == d && rend == s && (d ? ps + piece == nextref : ps == nextref);
                    if (!cont) { in = true; rdir = d; rstart = s; }
                    rend = s + piece; nextref = d ? ps : ps + piece;
                    const uint64_t lim = std::min(s + piece, rstart + K - 1); if (s < lim) lookup_range(s, lim);
                    if (!d) src += piece;
                } else { in = false; lookup_range(s, s + piece); }
                s += piece; len -= piece; } } }
}

static int cmd_residual(char** argv) {
    Ref ref = load_ref(argv[2]); const std::string dir = argv[3]; std::vector<std::string> names = read_names(argv[4]);
    Matcher G; G.M = atoi(argv[6]); std::vector<std::string> gn;
    { FILE* f = fopen(argv[5], "r"); char b[4096]; fgets(b, sizeof b, f); while (fgets(b, sizeof b, f)) { char so[64], sq[64]; unsigned gid;
        if (sscanf(b, "%u\t%63s\t%63s", &gid, so, sq) == 3) { uint64_t x = 0; for (int i = 0; i < 20; i++) x = (x << 2) | (uint64_t)b2((uint8_t)sq[i]); G.fw.push_back(x); G.rv.push_back(rc_n(x, 20)); gn.push_back(so); } } fclose(f); }
    const int threads = atoi(argv[7]); const double t0 = now_s();
    std::vector<ResMap> Rth(threads); std::vector<std::vector<std::string>> cname(names.size());
    std::atomic<size_t> next{0}; std::vector<std::thread> th;
    for (int t = 0; t < threads; t++) th.emplace_back([&, t] { for (;;) { const size_t i = next.fetch_add(1); if (i >= names.size()) break;
        const std::string p = dir + "/" + names[i] + ".q16k.rr3"; Archive X; open_v1(slurp(p), p, ref.sha, ref.n, X);
        std::set<std::pair<uint32_t, uint64_t>> seen; search_residual(X, ref, G, (uint32_t)i, Rth[t], seen);
        std::vector<std::string> cn; for (auto& c : X.rec) cn.push_back(c.name); cname[i] = std::move(cn); } });
    for (auto& x : th) x.join();
    ResMap R; for (auto& m : Rth) for (auto& kv : m) { ResRec& r = R[kv.first];
        if (!r.occ) { r = kv.second; } else { r.occ += kv.second.occ; r.asms += kv.second.asms; } }
    FILE* o = fopen(argv[8], "w"); fprintf(o, "guide\tcanonical_seq\toccurrences\tassemblies\tmismatches\tex_assembly\tex_contig\tex_pos0\tex_strand\n");
    for (auto& kv : R) { const ResRec& r = kv.second;
        fprintf(o, "%s\t%s\t%llu\t%llu\t%u\t%s\t%s\t%llu\t%c\n", gn[kv.first.first].c_str(), code_str(kv.first.second).c_str(), (unsigned long long)r.occ,
                (unsigned long long)r.asms, r.mm, names[r.a].c_str(), r.c < cname[r.a].size() ? cname[r.a][r.c].c_str() : "?", (unsigned long long)r.pos, r.strand ? '-' : '+'); }
    fclose(o);
    fprintf(stderr, "RESIDUAL\tassemblies %zu\tguides %zu\tunique sequences %zu\twall %.1f s\tthreads %d\tpeak RSS %.1f GB\n", names.size(), G.fw.size(), R.size(), now_s() - t0, threads, peak_rss_gb());
    return 0;
}

// kmerscan: how many times each given 23-mer occurs in the reference (either strand), one scan
static int cmd_kmerscan(char** argv) {
    Ref ref = load_ref(argv[2]);
    std::vector<uint64_t> pat; std::vector<std::string> txt;
    { FILE* f = fopen(argv[3], "r"); char b[256]; while (fgets(b, sizeof b, f)) { std::string q(b); while (!q.empty() && (q.back() == '\n' || q.back() == '\r')) q.pop_back();
        if (q.size() != (size_t)K) continue; uint64_t x = 0; bool ok = true; for (int i = 0; i < K; i++) { const int v = b2((uint8_t)q[i]); if (v < 0) { ok = false; break; } x = (x << 2) | (uint64_t)v; }
        if (ok) { pat.push_back(x); txt.push_back(q); } } fclose(f); }
    std::map<uint64_t, size_t> idx; for (size_t i = 0; i < pat.size(); i++) { idx[pat[i]] = i; idx[rc23(pat[i])] = i; }
    const int threads = atoi(argv[4]); std::vector<std::vector<uint64_t>> cnt(threads, std::vector<uint64_t>(pat.size(), 0));
    const uint64_t n = ref.n, chunk = (n + threads - 1) / threads; std::vector<std::thread> th; const double t0 = now_s();
    for (int t = 0; t < threads; t++) th.emplace_back([&, t] { const uint64_t a = t * chunk, b = std::min(n, a + chunk); Roller r;
        for (uint64_t p = a >= (uint64_t)(K - 1) ? a - (K - 1) : 0; p < std::min<uint64_t>(n, b + K - 1); p++) if (r.push(ref.R[p])) { const uint64_t st = p + 1 - K;
            if (st >= a && st < b) { auto it = idx.find(r.code); if (it != idx.end()) cnt[t][it->second]++; } } });
    for (auto& x : th) x.join();
    FILE* o = fopen(argv[5], "w"); fprintf(o, "seq\treference_occurrences\n");
    for (size_t i = 0; i < pat.size(); i++) { uint64_t c = 0; for (int t = 0; t < threads; t++) c += cnt[t][i]; fprintf(o, "%s\t%llu\n", txt[i].c_str(), (unsigned long long)c); }
    fclose(o);
    fprintf(stderr, "KMERSCAN\tsequences %zu\treference bases %llu\tscan %.1f s\tthreads %d\tpeak RSS %.1f GB\n", pat.size(), (unsigned long long)n, now_s() - t0, threads, peak_rss_gb());
    return 0;
}

// origin: for each guide, split its sites into those transferred from the single reference scan (the window lies
// inside one reference-contiguous run) and those found in the assembly's own residual windows (literals, self copies,
// run boundaries). For residual sites the context +/- 20 bases is hashed so unique sequences can be counted.
struct GOrigin { std::vector<uint64_t> total, transfer, residual, asms; std::vector<std::array<uint64_t, 8>> bymm;
    std::vector<std::vector<uint64_t>> ctx;
    void init(size_t n) { total.assign(n, 0); transfer.assign(n, 0); residual.assign(n, 0); asms.assign(n, 0); bymm.assign(n, {}); ctx.resize(n); } };
static const int GPRE = 64;

static void search_gorigin(const Archive& X, const Ref& ref, const Matcher& G, const std::vector<RefSite>& occ, uint64_t& nsites, GOrigin& O, std::vector<uint64_t>& res_local) {
    std::vector<RrOp> ops(RR_MAXOPS); std::vector<uint8_t> lit(X.Q), ext(GPRE + X.Q, 'N');
    std::vector<uint64_t> cstart, cend; for (auto& c : X.rec) { cstart.push_back(c.boff); cend.push_back(c.boff + c.len); }
    size_t ci = 0; while (ci < cend.size() && cend[ci] == 0) ci++;
    std::vector<Run> runs; bool in = false; int rdir = 0; uint64_t rstart = 0, rend = 0, nextref = 0, rlo = 0;
    auto close_run = [&]() { if (in && rend - rstart >= (uint64_t)K) runs.push_back({rstart, rend - rstart, rlo, rdir, (uint32_t)ci}); in = false; };
    uint64_t last_e = ~0ull, code = 0; int valid = 0; uint64_t bs = 0; uint32_t blen = 0;
    struct Pend { uint64_t e; uint32_t ci, gid; }; std::vector<Pend> pend, pend_next;
    auto hash_ctx = [&](uint64_t e, uint32_t c) -> uint64_t {
        const uint64_t p = e + 1 - K, lo = std::max<uint64_t>(cstart[c], p >= 20 ? p - 20 : 0);
        const uint64_t hi = std::min<uint64_t>(std::min<uint64_t>(cend[c], e + 21), bs + blen);
        return XXH3_64bits(ext.data() + GPRE + ((int64_t)lo - (int64_t)bs), (size_t)(hi - lo)); };
    for (uint64_t b = 0; b < X.nb; b++) {
        bs = b * X.Q; blen = (uint32_t)std::min<uint64_t>(X.Q, X.nbases - bs);
        uint32_t o, l; R3Diag g; X.locate(b, &o, &l, &g);
        const int k = r3_decode_block(X.P + o, l, X.T, ref.R.data(), ref.n, blen, g, ops.data(), (uint32_t)ops.size(), lit.data(), (uint32_t)lit.size());
        if (k < 0) throw Err(X.path + ": decode");
        if (b) memmove(ext.data(), ext.data() + X.Q, GPRE);
        exec_ops(ops.data(), k, lit.data(), ref.R.data(), ext.data() + GPRE);
        for (const Pend& q : pend) O.ctx[q.gid].push_back(hash_ctx(q.e, q.ci));
        pend.clear();
        const uint8_t* base = ext.data() + GPRE;
        auto lookup_range = [&](uint64_t e0, uint64_t e1) {
            for (uint64_t e = e0; e < e1; e++) {
                if (last_e == ~0ull || last_e + 1 != e) {
                    const uint64_t x0 = std::max<uint64_t>(cstart[ci], e + 1 >= (uint64_t)K ? e + 1 - K : 0); code = 0; valid = 0;
                    for (uint64_t x = x0; x <= e; x++) { const int v = b2(base[(int64_t)x - (int64_t)bs]); if (v < 0) { valid = 0; code = 0; } else { code = ((code << 2) | (uint64_t)v) & KMASK; valid++; } }
                } else { const int v = b2(base[(int64_t)e - (int64_t)bs]); if (v < 0) { valid = 0; code = 0; } else { code = ((code << 2) | (uint64_t)v) & KMASK; if (valid < K) valid++; } }
                last_e = e;
                if (valid >= K) G.match(code, [&](uint32_t gid, int, int mm) { nsites++; O.total[gid]++; O.residual[gid]++; O.bymm[gid][mm]++; res_local[gid]++;
                    if (e + 21 <= bs + blen) O.ctx[gid].push_back(hash_ctx(e, (uint32_t)ci)); else pend_next.push_back({e, (uint32_t)ci, gid}); }); } };
        for (int j = 0; j < k; j++) { const RrOp& q = ops[j]; uint64_t s = bs + q.dst, len = q.len; uint64_t src = q.src;
            while (len) {
                while (ci < cend.size() && s >= cend[ci]) { close_run(); ci++; last_e = ~0ull; }
                const uint64_t piece = std::min<uint64_t>(len, cend[ci] - s);
                if (q.kind == 1 || q.kind == 3) { const int d = q.kind == 3; const uint64_t ps = d ? src + len - piece : src;
                    const bool cont = in && rdir == d && rend == s && (d ? ps + piece == nextref : ps == nextref);
                    if (!cont) { close_run(); in = true; rdir = d; rstart = s; rlo = ps; }
                    rend = s + piece; nextref = d ? ps : ps + piece; if (d) rlo = ps;
                    const uint64_t lim = std::min(s + piece, rstart + K - 1); if (s < lim) lookup_range(s, lim);
                    if (!d) src += piece;
                } else { close_run(); lookup_range(s, s + piece); }
                s += piece; len -= piece; } }
        pend.swap(pend_next); pend_next.clear(); }
    for (const Pend& q : pend) O.ctx[q.gid].push_back(hash_ctx(q.e, q.ci));
    while (ci < cend.size()) { close_run(); ci++; }
    std::sort(runs.begin(), runs.end(), [](const Run& x, const Run& y) { return x.ref_lo < y.ref_lo; });
    size_t lo = 0;
    for (const Run& R_ : runs) { while (lo < occ.size() && occ[lo].pos < R_.ref_lo) lo++; const uint64_t last = R_.ref_lo + R_.len - K;
        for (size_t j = lo; j < occ.size() && occ[j].pos <= last; j++) { nsites++; O.total[occ[j].gid]++; O.transfer[occ[j].gid]++; O.bymm[occ[j].gid][occ[j].mm]++; } }
}

static int cmd_gorigin(char** argv) {
    Ref ref = load_ref(argv[2]); const std::string dir = argv[3]; std::vector<std::string> names = read_names(argv[4]);
    Matcher G; G.M = atoi(argv[6]); std::vector<std::string> gn;
    { FILE* f = fopen(argv[5], "r"); char b[4096]; fgets(b, sizeof b, f); while (fgets(b, sizeof b, f)) { char so[64], sq[64]; unsigned gid;
        if (sscanf(b, "%u\t%63s\t%63s", &gid, so, sq) == 3) { uint64_t x = 0; for (int i = 0; i < 20; i++) x = (x << 2) | (uint64_t)b2((uint8_t)sq[i]); G.fw.push_back(x); G.rv.push_back(rc_n(x, 20)); gn.push_back(so); } } fclose(f); }
    const int threads = atoi(argv[7]); const size_t NG = G.fw.size(); const double t0 = now_s();
    std::vector<RefSite> occ = scan_reference(ref, G, threads); const double tref = now_s() - t0;
    std::vector<GOrigin> Oth(threads); for (auto& o : Oth) o.init(NG);
    std::vector<uint64_t> ns(threads, 0); std::atomic<size_t> next{0}; std::vector<std::thread> th; const double t1 = now_s();
    for (int t = 0; t < threads; t++) th.emplace_back([&, t] { std::vector<uint64_t> rl(NG, 0);
        for (;;) { const size_t i = next.fetch_add(1); if (i >= names.size()) break;
            const std::string p = dir + "/" + names[i] + ".q16k.rr3"; Archive X; open_v1(slurp(p), p, ref.sha, ref.n, X);
            std::fill(rl.begin(), rl.end(), 0); search_gorigin(X, ref, G, occ, ns[t], Oth[t], rl);
            for (size_t q = 0; q < NG; q++) if (rl[q]) Oth[t].asms[q]++; } });
    for (auto& x : th) x.join(); const double tasm = now_s() - t1;
    GOrigin O; O.init(NG); uint64_t nsites = 0;
    for (auto& o : Oth) { nsites += 0; for (size_t i = 0; i < NG; i++) { O.total[i] += o.total[i]; O.transfer[i] += o.transfer[i]; O.residual[i] += o.residual[i]; O.asms[i] += o.asms[i]; for (int m = 0; m < 8; m++) O.bymm[i][m] += o.bymm[i][m]; } }
    for (auto& v : ns) nsites += v;
    printf("GORIGIN\tassemblies %zu\tguides %zu\tmax mismatches %d\tsites %llu\treference scan %.1f s\tassemblies %.1f s\ttotal %.1f s\tthreads %d\tpeak RSS %.1f GB\n",
           names.size(), NG, G.M, (unsigned long long)nsites, tref, tasm, now_s() - t0, threads, peak_rss_gb());
    printf("guide\ttotal\ttransfer\tresidual\tresidual_frac\tassemblies_with_residual\tunique_contexts\tmm0\tmm1\tmm2\tmm3\tmm0_3\n");
    for (size_t i = 0; i < NG; i++) {
        std::vector<uint64_t> buf; for (auto& o : Oth) { buf.insert(buf.end(), o.ctx[i].begin(), o.ctx[i].end()); std::vector<uint64_t>().swap(o.ctx[i]); }
        std::sort(buf.begin(), buf.end()); const uint64_t u = (uint64_t)(std::unique(buf.begin(), buf.end()) - buf.begin());
        const uint64_t m03 = O.bymm[i][0] + O.bymm[i][1] + O.bymm[i][2] + O.bymm[i][3];
        printf("%s\t%llu\t%llu\t%llu\t%.6f\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\n", gn[i].c_str(), (unsigned long long)O.total[i], (unsigned long long)O.transfer[i],
               (unsigned long long)O.residual[i], O.total[i] ? (double)O.residual[i] / (double)O.total[i] : 0.0, (unsigned long long)O.asms[i], (unsigned long long)u,
               (unsigned long long)O.bymm[i][0], (unsigned long long)O.bymm[i][1], (unsigned long long)O.bymm[i][2], (unsigned long long)O.bymm[i][3], (unsigned long long)m03); }
    return 0;
}

// sites: dump every site of every assembly (parse path), with contig names, for comparison with another tool
static int cmd_sites(char** argv) {
    Ref ref = load_ref(argv[2]); const std::string dir = argv[3]; std::vector<std::string> names = read_names(argv[4]);
    Matcher G; G.M = atoi(argv[6]); std::vector<std::string> gn;
    { FILE* f = fopen(argv[5], "r"); char b[4096]; fgets(b, sizeof b, f); while (fgets(b, sizeof b, f)) { char so[64], sq[64]; unsigned gid;
        if (sscanf(b, "%u\t%63s\t%63s", &gid, so, sq) == 3) { uint64_t x = 0; for (int i = 0; i < 20; i++) x = (x << 2) | (uint64_t)b2((uint8_t)sq[i]); G.fw.push_back(x); G.rv.push_back(rc_n(x, 20)); gn.push_back(so); } } fclose(f); }
    const int threads = atoi(argv[7]); const double t0 = now_s();
    std::vector<RefSite> occ = scan_reference(ref, G, threads); const double tref = now_s() - t0;
    const size_t N = names.size(); std::vector<std::string> outs(N); std::atomic<size_t> next{0}; std::vector<Stage> S(threads); std::vector<std::thread> th;
    const double t1 = now_s();
    for (int t = 0; t < threads; t++) th.emplace_back([&, t] { for (;;) { const size_t i = next.fetch_add(1); if (i >= N) break;
        const std::string p = dir + "/" + names[i] + ".q16k.rr3"; Archive X; open_v1(slurp(p), p, ref.sha, ref.n, X);
        std::vector<Hit> H; search_parse(X, ref, G, occ, H, S[t]); std::sort(H.begin(), H.end());
        std::string o; char buf[512];
        for (const Hit& h : H) { snprintf(buf, sizeof buf, "%s\t%s\t%llu\t%c\t%u\t%s\n", names[i].c_str(), X.rec[h.contig].name.c_str(),
            (unsigned long long)h.pos, h.strand ? '-' : '+', h.mm, gn[h.gid].c_str()); o += buf; }
        outs[i] = std::move(o); } });
    for (auto& x : th) x.join(); const double tasm = now_s() - t1;
    FILE* o = fopen(argv[8], "w"); fprintf(o, "assembly\tcontig\tpos0\tstrand\tmismatches\tguide\n"); uint64_t n = 0;
    for (auto& v : outs) { fwrite(v.data(), 1, v.size(), o); for (char c : v) n += c == '\n'; } fclose(o);
    fprintf(stderr, "SITES\tassemblies %zu\tguides %zu\tmax mismatches %d\tsites %llu\treference scan %.1f s\tassemblies %.1f s\ttotal %.1f s\tthreads %d\tpeak RSS %.1f GB\n",
            N, G.fw.size(), G.M, (unsigned long long)n, tref, tasm, now_s() - t0, threads, peak_rss_gb());
    return 0;
}

// refscan: dump the reference sites of every guide (one scan of the reference), for comparison with another tool
static int cmd_refscan(char** argv) {
    Ref ref = load_ref(argv[2]); Matcher G; G.M = atoi(argv[4]); std::vector<std::string> nm;
    { FILE* f = fopen(argv[3], "r"); char b[4096]; fgets(b, sizeof b, f); while (fgets(b, sizeof b, f)) { char so[64], sq[64]; unsigned gid;
        if (sscanf(b, "%u\t%63s\t%63s", &gid, so, sq) == 3) { uint64_t x = 0; for (int i = 0; i < 20; i++) x = (x << 2) | (uint64_t)b2((uint8_t)sq[i]); G.fw.push_back(x); G.rv.push_back(rc_n(x, 20)); nm.push_back(so); } } fclose(f); }
    const int threads = atoi(argv[5]); const double t0 = now_s();
    std::vector<RefSite> occ = scan_reference(ref, G, threads); const double dt = now_s() - t0;
    FILE* o = fopen(argv[6], "w"); fprintf(o, "guide\tpos0\tstrand\tmismatches\tsite\n");
    for (const RefSite& x : occ) { std::string w((const char*)&ref.R[x.pos], K); fprintf(o, "%s\t%llu\t%c\t%u\t%s\n", nm[x.gid].c_str(), (unsigned long long)x.pos, x.strand ? '-' : '+', x.mm, w.c_str()); }
    fclose(o);
    fprintf(stderr, "REFSCAN\tguides %zu\tmax mismatches %d\treference bases %llu\tsites %zu\tscan %.1f s\tthreads %d\tpeak RSS %.1f GB\n",
            G.fw.size(), G.M, (unsigned long long)ref.n, occ.size(), dt, threads, peak_rss_gb());
    return 0;
}

int main(int argc, char** argv) {
    try { if (argc >= 9 && !strcmp(argv[1], "residual")) return cmd_residual(argv);
          if (argc >= 6 && !strcmp(argv[1], "kmerscan")) return cmd_kmerscan(argv);
          if (argc >= 8 && !strcmp(argv[1], "gorigin")) return cmd_gorigin(argv);
          if (argc >= 9 && !strcmp(argv[1], "sites")) return cmd_sites(argv);
          if (argc >= 7 && !strcmp(argv[1], "refscan")) return cmd_refscan(argv);
          if (argc >= 10 && !strcmp(argv[1], "guides")) return cmd_guides(argv); if (argc >= 9 && !strcmp(argv[1], "run")) return cmd_run(argc, argv); }
    catch (const std::exception& e) { fprintf(stderr, "%s\n", e.what()); return 2; }
    fprintf(stderr, "usage: see header\n"); return 1;
}
