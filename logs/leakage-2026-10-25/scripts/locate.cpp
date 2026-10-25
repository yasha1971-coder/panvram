// rv2dec.cpp - refrel v2 PROTOTYPE, DECODER ONLY (H4 cohort anchors; format NOT frozen, marked unstable; v1 untouched).
// This is research/refrel_v2/rv2.cpp of the aceapex refrel-v2 branch with the encoder removed (the v2 parse, symbol
// generation, rANS encoding, block-table writing and the commands that call them); every decode path is unchanged.
// Grammar: rv2.h. Container "RFRL3P2": the v1 header layout (FORMAT.md section 1) with magic RFRL3P2\0 and version 2
// or 3, then the reference name, then the ANCHOR TABLE (u16 K; per anchor u16 name length, name, u64 bases, u64 XXH3 of
// its decoded upper-case bases), then the meta (zstd; v1 layout with 66 model contexts; version 3 = bt3 block table),
// the block-hash section (flag bit 0; not used here), the payload. Header XXH3 covers [0, payload start) as in v1.
// Anchors themselves are refrel3 v1 archives (depth 1).
//   judge  <t2t.fa> <v1 dir> <anchors|-> <rv2 dir> <manifest_v1.tsv> <names.txt> <threads>   every name decoded from
//          its archive (anchor: the v1 archive; else the rv2 archive) to FASTA and compared with the source FASTA XXH3
//          and size of the manifest; sizes (anchors: v1 bytes without the block-hash section)
//   bench  req|each <t2t.fa> <v1 dir> <anchors|-> <rv2 dir> <requests.tsv|samples> <truth.tsv> <threads> <warmups> <runs> <label>
//          as rr3bench (anchors from their v1 archives with T2T only, others rv2)
//   info   <in.rv2>
//   ledger, absstudy, wintruth, lawtruth, law: see each command below
// Build: research/refrel_v2/build.sh
#include "refrel3v1_nomain.cpp"
#include "rv2.h"
#include "bench_common.h"
#include <memory>
#include <sys/stat.h>

static const char P2_MAGIC[8] = {'R', 'F', 'R', 'L', '3', 'P', '2', 0};
static const uint32_t P2_VERSION = 2;          // block table interleaved (length, state = zigzag(c - (prev c +- Q)) x 2 + dir)
static const uint32_t P3_VERSION = 3;          // lever 1: block table columnar with per-target state prediction (bt3 below)
// bt3: per block, in three columns: LEB payload length; one target-code byte (0 = the target of the previous block's
// start state, 1 = outside every target, 2 + t = target t); LEB zigzag(x - pred) x 2 + dir, x = the start diagonal c,
// pred = the last start state of that target advanced by Q per block in its strand (the target base if none, 0 outside)
static int bt3_target(const std::vector<uint64_t>& base, const std::vector<uint64_t>& len, int64_t x) {
    const int NT = (int)base.size(); if (x < 0) return NT; for (int t = NT - 1; t >= 0; t--) if ((uint64_t)x >= base[t]) return (uint64_t)x < base[t] + len[t] ? t : NT; return NT; }

// ---------------------------------------------------------------- cohort: T2T + anchors (decoded from v1) + lift tables
struct Cohort { Ref ref; std::vector<uint8_t> G; Space S; std::vector<std::string> an; std::vector<uint64_t> axxh; };
static std::vector<std::string> split_c(const std::string& s) { std::vector<std::string> v; if (s == "-" || s.empty()) return v; size_t p = 0; for (;;) { size_t c = s.find(',', p); v.push_back(s.substr(p, c - p)); if (c == std::string::npos) break; p = c + 1; } return v; }
static std::string v1path(const std::string& dir, const std::string& nm) { return dir + "/" + nm + ".q16k.rr3"; }
// T2T + every anchor decoded from its v1 archive (all v1 open checks, block XXH3 checked) into one buffer; lift tables
// = each anchor's reference copies in stream order (+ index sorted by (p, a))
static void load_cohort(Cohort& C, const char* t2t, const std::string& v1dir, const std::vector<std::string>& an, int T) {
    C.ref = load_ref(t2t, 1, false); C.an = an; const uint64_t Rn = C.ref.R.size();
    std::vector<std::vector<uint8_t>> files(an.size()); uint64_t total = Rn;
    for (size_t k = 0; k < an.size(); k++) { files[k] = slurp(v1path(v1dir, an[k])); total += rd64(&files[k][64]); }
    C.G.reserve(total); C.G.assign(C.ref.R.begin(), C.ref.R.end()); std::vector<uint8_t>().swap(C.ref.R);
    C.S.t.resize(1 + an.size()); C.S.t[0].base = 0; C.S.t[0].len = Rn;
    Ref view; memcpy(view.sha, C.ref.sha, 32); view.name = C.ref.name; if (!an.empty()) view.R.assign(C.G.begin(), C.G.begin() + Rn);   // open_v1 checks size + SHA against a Ref
    for (size_t k = 0; k < an.size(); k++) {
        V1 X; std::string why; if (open_v1(files[k].data(), files[k].size(), view, X, why)) { fprintf(stderr, "anchor %s refused: %s\n", an[k].c_str(), why.c_str()); exit(2); }
        const uint64_t base = C.G.size(); C.G.resize(base + X.nbases); Tgt& tg = C.S.t[1 + k]; tg.base = base; tg.len = X.nbases;
        std::vector<std::vector<LiftE>> le(X.nb); std::atomic<uint64_t> next{0}; std::atomic<int> err{0}; std::vector<std::thread> th;
        for (int t = 0; t < T; t++) th.emplace_back([&] { std::vector<RrOp> ops(RR_MAXOPS); std::vector<uint8_t> lit(X.Q);
            for (;;) { const uint64_t b = next.fetch_add(1); if (b >= X.nb) break; const uint32_t blen = (uint32_t)std::min<uint64_t>(X.Q, X.nbases - b * X.Q);
                const int n = r3_decode_block(X.P + X.off[b], (uint32_t)(X.off[b + 1] - X.off[b]), X.T, C.G.data(), Rn, blen, X.st[b], ops.data(), (uint32_t)ops.size(), lit.data(), (uint32_t)lit.size());
                if (n < 0) { err++; continue; } uint8_t* out = &C.G[base + b * X.Q]; exec_fast(ops.data(), n, lit.data(), C.G.data(), out);
                if (X.hashes && XXH3_64bits(out, blen) != rd64(X.hashes + 8 * b)) err++;
                for (int j = 0; j < n; j++) if (ops[j].kind == 1 || ops[j].kind == 3) le[b].push_back({ops[j].src, b * X.Q + ops[j].dst, ops[j].len, ops[j].kind == 3 ? 1u : 0u}); } });
        for (auto& x : th) x.join(); if (err) { fprintf(stderr, "anchor %s: %d blocks failed\n", an[k].c_str(), (int)err); exit(3); }
        size_t ne = 0; for (auto& v : le) ne += v.size(); tg.byA.reserve(ne); for (auto& v : le) { tg.byA.insert(tg.byA.end(), v.begin(), v.end()); std::vector<LiftE>().swap(v); }
        tg.byP.resize(ne); for (size_t i = 0; i < ne; i++) tg.byP[i] = (uint32_t)i;
        std::sort(tg.byP.begin(), tg.byP.end(), [&](uint32_t x, uint32_t y) { const LiftE &a = tg.byA[x], &b = tg.byA[y]; return a.p != b.p ? a.p < b.p : a.a < b.a; });
        C.axxh.push_back(XXH3_64bits(&C.G[base], X.nbases));
        fprintf(stderr, "anchor %s: %llu bases, %zu lift entries\n", an[k].c_str(), (unsigned long long)X.nbases, ne);
    }
    std::vector<uint8_t>().swap(view.R);
    C.S.G = C.G.data(); C.S.Gn = C.G.size();
}

static void build_tab4(R4Tab& T) { for (int c = 0; c < R4_NCTX; c++) { uint32_t a = 0; for (int s = 0; s < R3_AB; s++) { T.cum[c][s] = (uint16_t)a; for (uint32_t k = 0; k < T.freq[c][s]; k++) T.sym[c][a + k] = (uint8_t)s; a += T.freq[c][s]; } T.cum[c][R3_AB] = (uint16_t)a; } }

// ---------------------------------------------------------------- reader: every check of v1 + the anchor table against the loaded cohort
struct V2 { uint32_t Q = 0; uint64_t nbases = 0, nb = 0, fasta_xxh = 0; std::vector<Rec> rec; std::vector<std::pair<uint64_t, uint64_t>> low; R4Tab* T = nullptr;
    std::vector<uint64_t> off; std::vector<R3Diag> st; const uint8_t* P = nullptr; ~V2() { delete T; } };
static int open_v2(const uint8_t* a, size_t n, const Cohort& C, V2& X, std::string& why) {
    auto fail = [&](const char* w) { why = w; return 1; };
    if (n < V1_HDR + 4) return fail("short file");
    if (memcmp(a, P2_MAGIC, 8)) return fail("magic");
    const uint32_t ver = rd32(a + 8); if (ver != P2_VERSION && ver != P3_VERSION) return fail("version");
    X.Q = rd32(a + 12); if (!q_ok(X.Q)) return fail("block size");
    if (rd32(a + 16) || rd32(a + 20)) return fail("flags / reserved");
    X.nbases = rd64(a + 64); X.nb = rd64(a + 72); X.fasta_xxh = rd64(a + 80);
    const uint64_t mz = rd64(a + 96), mr = rd64(a + 104), hs = rd64(a + 112), pl = rd64(a + 120);
    if (hs) return fail("hash section");
    size_t q = V1_HDR; const uint16_t nl = (uint16_t)(a[q] | a[q + 1] << 8); q += 2 + nl; if (q + 2 > n) return fail("short file");
    const uint16_t K_ = (uint16_t)(a[q] | a[q + 1] << 8); q += 2;
    if (K_ != C.an.size()) return fail("anchor count differs from the loaded cohort");
    for (size_t k = 0; k < K_; k++) { if (q + 2 > n) return fail("anchor table"); const uint16_t l = (uint16_t)(a[q] | a[q + 1] << 8); q += 2; if (q + l + 16 > n) return fail("anchor table");
        if (std::string((const char*)a + q, l) != C.an[k]) return fail("anchor name differs"); q += l;
        if (rd64(a + q) != C.S.t[1 + k].len || rd64(a + q + 8) != C.axxh[k]) return fail("anchor bases differ"); q += 16; }
    const uint64_t pstart = q + mz;
    if (mz > n || pl > n || pstart > n || pstart + pl != n) return fail("section sizes");
    if (X.nb != (X.nbases + X.Q - 1) / X.Q || (X.nbases && !X.nb)) return fail("block count");
    { std::vector<uint8_t> h(a, a + pstart); memset(&h[128], 0, 8); if (XXH3_64bits(h.data(), h.size()) != rd64(a + 128)) return fail("header XXH3"); }
    if (memcmp(a + 24, C.ref.sha, 32) || rd64(a + 56) != C.S.t[0].len) return fail("reference SHA-256 / size differs");
    const uint8_t* mzp = a + q;
    if (mr > (1ull << 32) || ZSTD_getFrameContentSize(mzp, mz) != mr) return fail("meta frame");
    std::vector<uint8_t> M(mr); if (ZSTD_decompress(M.data(), mr, mzp, mz) != mr) return fail("meta decompress");
    size_t i = 0; bool bad = false;
    auto need = [&](size_t k) { if (i + k > M.size()) { bad = true; return false; } return true; };
    auto g64 = [&]() -> uint64_t { if (!need(8)) return 0; uint64_t v; memcpy(&v, &M[i], 8); i += 8; return v; };
    auto g32 = [&]() -> uint32_t { if (!need(4)) return 0; uint32_t v; memcpy(&v, &M[i], 4); i += 4; return v; };
    auto gl = [&]() -> uint64_t { uint64_t v = 0; int sh = 0; for (;;) { if (!need(1) || sh > 63) { bad = true; return 0; } uint8_t c = M[i++]; v |= (uint64_t)(c & 0x7F) << sh; if (!(c & 0x80)) break; sh += 7; } return v; };
    const uint32_t nr = g32(); uint64_t bo = 0; if (nr > M.size()) return fail("records");
    for (uint32_t r = 0; r < nr && !bad; r++) { Rec rr; const uint32_t hl = g32(); if (!need(hl)) break; rr.hdr.assign((const char*)&M[i], hl); i += hl; rr.len = g64(); rr.lw = g32(); if (!rr.lw && rr.len) bad = true; rr.boff = bo; bo += rr.len; rr.foff = 0; X.rec.push_back(rr); }
    if (bad || bo != X.nbases) return fail("contig table");
    const uint64_t nlw = g64(); uint64_t lst = 0; if (nlw > X.nbases) return fail("case runs");
    for (uint64_t k = 0; k < nlw && !bad; k++) { const uint64_t g = gl(), l = gl(); if (lst + g + l > X.nbases) { bad = true; break; } X.low.push_back({lst + g, l}); lst += g + l; }
    if (bad) return fail("case runs");
    X.T = new R4Tab(); memset(X.T, 0, sizeof(R4Tab));
    for (int c = 0; c < R4_NCTX && !bad; c++) { uint64_t sum = 0; for (int s = 0; s < R4_ALPHA[c]; s++) { const uint64_t f = gl(); if (f > R3_M) bad = true; X.T->freq[c][s] = (uint16_t)f; sum += f; } if (sum && sum != R3_M) bad = true; }
    if (bad) return fail("model tables");
    build_tab4(*X.T);
    X.off.assign(X.nb + 1, 0); X.st.resize(X.nb); uint64_t prevc = 0; uint32_t prevd = 0;
    if (ver == P2_VERSION) {
    for (uint64_t b = 0; b < X.nb && !bad; b++) { const uint64_t len = gl(); X.off[b + 1] = X.off[b] + len; const uint64_t v = gl();
        const int64_t pred = b == 0 ? 0 : (int64_t)(prevd ? prevc - X.Q : prevc + X.Q); const uint64_t qq = v >> 1;
        X.st[b].dir = (uint32_t)(v & 1); X.st[b].c = (uint64_t)(pred + ((int64_t)(qq >> 1) ^ -(int64_t)(qq & 1))); prevc = X.st[b].c; prevd = X.st[b].dir; } }
    else {                                                                    // bt3 (lever 1)
        std::vector<uint64_t> base, len; for (auto& tg : C.S.t) { base.push_back(tg.base); len.push_back(tg.len); } const int NT = (int)base.size();
        for (uint64_t b = 0; b < X.nb && !bad; b++) { const uint64_t l = gl(); X.off[b + 1] = X.off[b] + l; }
        if (bad || i + X.nb > M.size()) return fail("block table");
        const uint8_t* tcp = &M[i]; i += X.nb;
        std::vector<int64_t> lx(NT + 1, 0); std::vector<uint64_t> lb(NT + 1, 0); std::vector<uint32_t> ld(NT + 1, 0); std::vector<uint8_t> seen(NT + 1, 0); int tprev = -1;
        for (uint64_t b = 0; b < X.nb && !bad; b++) { const uint8_t code = tcp[b]; int t;
            if (code == 0) { if (tprev < 0) return fail("block table"); t = tprev; } else if (code == 1) t = NT; else { t = code - 2; if (t >= NT) return fail("block table"); }
            const uint64_t v = gl(); const uint64_t qq = v >> 1;
            const int64_t pred = (seen[t] ? lx[t] + (int64_t)(b - lb[t]) * (ld[t] ? -(int64_t)X.Q : (int64_t)X.Q) : (t == NT ? 0 : (int64_t)base[t]))
#ifdef RV2_FAULT_BT3
                + 1                                                 // negative control: bt3 per-target prediction one base off
#endif
                ;
            const int64_t x = pred + ((int64_t)(qq >> 1) ^ -(int64_t)(qq & 1)); if (bt3_target(base, len, x) != t) return fail("block table target");
            X.st[b].dir = (uint32_t)(v & 1); X.st[b].c = (uint64_t)x; lx[t] = x; lb[t] = b; ld[t] = X.st[b].dir; seen[t] = 1; tprev = t; } }
    if (bad || i != M.size() || X.off[X.nb] != pl) return fail("block table");
    X.P = a + pstart; return 0;
}
static inline int block_v2(const V2& X, const Cohort& C, uint64_t b, uint8_t* out, std::vector<RrOp>& ops, std::vector<uint8_t>& lit) {
    const uint32_t blen = (uint32_t)std::min<uint64_t>(X.Q, X.nbases - b * X.Q);
    const int k = r4_decode_block(X.P + X.off[b], (uint32_t)(X.off[b + 1] - X.off[b]), X.T, C.S, blen, X.st[b], ops.data(), (uint32_t)ops.size(), lit.data(), (uint32_t)lit.size());
    if (k < 0) return 1;
    exec_fast(ops.data(), k, lit.data(), C.S.G, out); return 0;
}
static int full_v2(const V2& X, const Cohort& C, std::string& fa, int T) {
    std::vector<uint8_t> A(X.nbases + 64); std::atomic<uint64_t> next{0}; std::atomic<int> err{0};
    std::vector<std::thread> th; for (int t = 0; t < T; t++) th.emplace_back([&] { std::vector<RrOp> ops(RR_MAXOPS); std::vector<uint8_t> lit(X.Q);
        for (;;) { const uint64_t b = next.fetch_add(1); if (b >= X.nb || err) break; if (block_v2(X, C, b, &A[b * X.Q], ops, lit)) err = 1; } });
    for (auto& x : th) x.join(); if (err) return 1;
    for (auto& r : X.low) for (uint64_t x = r.first; x < r.first + r.second; x++) A[x] |= 0x20;
    fa.clear(); fa.reserve(X.nbases + X.nbases / 60 + 1024);
    for (auto& r : X.rec) { fa += '>'; fa += r.hdr; fa += '\n'; for (uint64_t x = 0; x < r.len; x += r.lw) { fa.append((const char*)&A[r.boff + x], std::min<uint64_t>(r.lw, r.len - x)); fa += '\n'; } }
    return XXH3_64bits(fa.data(), fa.size()) != X.fasta_xxh ? 3 : 0;
}

// ---------------------------------------------------------------- commands
static std::vector<std::string> read_names(const char* p) { std::vector<std::string> v; for (auto& r : read_tsv(p)) if (!r.empty() && !r[0].empty()) v.push_back(r[0]); return v; }
static int rv2_cmd_judge(char** argv) {
    const std::string v1dir = argv[3], dir = argv[5]; const std::vector<std::string> an = split_c(argv[4]); const int T = atoi(argv[8]);
    std::map<std::string, std::pair<uint64_t, uint64_t>> man;                               // name -> FASTA bytes, FASTA XXH3
    { auto rows = read_tsv(argv[6]); for (size_t r = 1; r < rows.size(); r++) man[rows[r][0]] = {strtoull(rows[r][4].c_str(), 0, 10), strtoull(rows[r][5].c_str(), 0, 16)}; }
    Cohort C; load_cohort(C, argv[2], v1dir, an, T);
    Ref view; memcpy(view.sha, C.ref.sha, 32); view.name = C.ref.name; view.R.assign(C.G.begin(), C.G.begin() + C.S.t[0].len);
    // optional (env RV2_TRUTH = truth.tsv with SAMPLE <set name> <bases> <SHA-256 of the upper-case bases>, env AVR_NAMES =
    // manifest name <tab> set name): the SHA-256 of the decoded upper-case bases must equal the truth as well
    std::map<std::string, std::string> tsha, set_of; if (getenv("RV2_TRUTH")) { for (auto& r : read_tsv(getenv("RV2_TRUTH"))) if (r[0] == "SAMPLE") tsha[r[1]] = r[3]; }
    if (getenv("AVR_NAMES")) for (auto& r : read_tsv(getenv("AVR_NAMES"))) set_of[r[0]] = r[1];
    auto sha_fa = [](const std::string& fa) { EVP_MD_CTX* c = EVP_MD_CTX_new(); EVP_DigestInit_ex(c, EVP_sha256(), nullptr); std::vector<uint8_t> buf; buf.reserve(1 << 20);
        size_t i = 0; while (i < fa.size()) { if (fa[i] == '>') { while (i < fa.size() && fa[i] != '\n') i++; i++; continue; } const char ch = fa[i++]; if (ch == '\n') continue;
            buf.push_back((uint8_t)((ch >= 'a' && ch <= 'z') ? ch - 32 : ch)); if (buf.size() == (1u << 20)) { EVP_DigestUpdate(c, buf.data(), buf.size()); buf.clear(); } }
        EVP_DigestUpdate(c, buf.data(), buf.size()); unsigned char md[32]; unsigned int l = 0; EVP_DigestFinal_ex(c, md, &l); EVP_MD_CTX_free(c); return sha_hex(md); };
    printf("JUDGE\tname\trole\tarchive_bytes\tcounted_bytes\tfasta_bytes\tfasta_xxh3\tresult\tdecode_s\tbases_sha256\n");
    int pass = 0, fail = 0; uint64_t tot = 0;
    for (auto& nm : read_names(argv[7])) { const bool isa = std::find(an.begin(), an.end(), nm) != an.end(); std::string fa, why; int e; uint64_t ab, cb; const double t0 = now_s();
        if (isa) { const std::vector<uint8_t> f = slurp(v1path(v1dir, nm)); V1 X; ab = f.size(); cb = f.size() - rd64(&f[112]);
            e = open_v1(f.data(), f.size(), view, X, why) ? 9 : full_v1(X, view, fa, true, T); }
        else { const std::vector<uint8_t> f = slurp(dir + "/" + nm + ".rv2"); V2 X; ab = cb = f.size();
            e = open_v2(f.data(), f.size(), C, X, why) ? 9 : full_v2(X, C, fa, T); }
        const double ds = now_s() - t0; const uint64_t x = XXH3_64bits(fa.data(), fa.size());
        bool ok = !e && man.count(nm) && man[nm].first == fa.size() && man[nm].second == x; std::string shs = "-";
        if (!tsha.empty()) { const std::string sn = set_of.count(nm) ? set_of[nm] : nm; const std::string h = e ? "" : sha_fa(fa); const bool sok = tsha.count(sn) && tsha[sn] == h;
            shs = (sok ? "sha256 ok " : "sha256 DIFFERS ") + h; ok = ok && sok; }
        ok ? pass++ : fail++; tot += cb;
        printf("JUDGE\t%s\t%s\t%llu\t%llu\t%zu\t%016llx\t%s\t%.1f\t%s\n", nm.c_str(), isa ? "anchor(v1)" : "rv2", (unsigned long long)ab, (unsigned long long)cb, fa.size(), (unsigned long long)x,
               ok ? "PASS" : (e == 9 ? ("REFUSED " + why).c_str() : e ? "DECODE-FAIL" : "DIFFERS"), ds, shs.c_str()); fflush(stdout); }
    printf("JUDGE_TOTAL\tpass %d\tfail %d\tcounted %llu B\t%.4f MB per sample\n", pass, fail, (unsigned long long)tot, tot / 1e6 / std::max(1, pass + fail));
    return fail ? 4 : 0;
}
static int rv2_cmd_info(char** argv) {
    const std::vector<uint8_t> a = slurp(argv[2]); if (a.size() < V1_HDR + 4 || memcmp(a.data(), P2_MAGIC, 8)) { printf("not rv2 prototype\n"); return 2; }
    size_t q = V1_HDR; const uint16_t nl = (uint16_t)(a[q] | a[q + 1] << 8); q += 2 + nl; const uint16_t K_ = (uint16_t)(a[q] | a[q + 1] << 8); q += 2;
    printf("rv2 PROTOTYPE v%u Q %u | reference %s | assembly %llu bases, %llu blocks | anchors %u:", rd32(&a[8]), rd32(&a[12]), std::string((const char*)&a[V1_HDR + 2], nl).c_str(),
           (unsigned long long)rd64(&a[64]), (unsigned long long)rd64(&a[72]), K_);
    for (int k = 0; k < K_; k++) { const uint16_t l = (uint16_t)(a[q] | a[q + 1] << 8); q += 2; printf(" %s", std::string((const char*)&a[q], l).c_str()); q += l + 16; }
    printf(" | meta %llu B, payload %llu B\n", (unsigned long long)rd64(&a[96]), (unsigned long long)rd64(&a[120])); return 0;
}
// bench: as rr3bench (research/agc_vs_refrel3/rr3bench.cpp); anchors from their v1 archives (T2T), others rv2
static int rv2_cmd_bench(int argc, char** argv) {
    if (argc < 13) { fprintf(stderr, "usage: see header\n"); return 1; }
    const std::string mode = argv[2], v1dir = argv[4], dir = argv[6]; const std::vector<std::string> an = split_c(argv[5]); Truth tr = read_truth(argv[8]);
    const int T = atoi(argv[9]), WU = atoi(argv[10]), R = atoi(argv[11]); const char* label = argv[12];
    const double t0 = now_s_(); Cohort C; load_cohort(C, argv[3], v1dir, an, 16);
    Ref view; memcpy(view.sha, C.ref.sha, 32); view.name = C.ref.name;                       // v1 path of the anchors: decode against G[0, T2T) (same bytes)
    std::vector<std::string> names; for (auto& kv : tr.bases) names.push_back(kv.first);
    std::map<std::string, size_t> idx; std::vector<std::vector<uint8_t>> bytes(names.size()); std::vector<std::unique_ptr<V1>> X1(names.size()); std::vector<std::unique_ptr<V2>> X2(names.size());
    std::vector<std::map<std::string, uint64_t>> boff(names.size()); std::vector<uint32_t> Qs(names.size());
    { view.R.assign(C.G.begin(), C.G.begin() + C.S.t[0].len);
      std::map<std::string, std::string> man; if (getenv("AVR_NAMES")) for (auto& v : read_tsv(getenv("AVR_NAMES"))) man[v[1]] = v[0];
      for (size_t i = 0; i < names.size(); i++) { idx[names[i]] = i; const std::string nm = man.count(names[i]) ? man[names[i]] : "y1_" + names[i]; std::string why; const bool isa = std::find(an.begin(), an.end(), nm) != an.end();
        const std::vector<Rec>* rec;
        if (isa) { bytes[i] = slurp(v1path(v1dir, nm)); X1[i].reset(new V1()); if (open_v1(bytes[i].data(), bytes[i].size(), view, *X1[i], why)) { fprintf(stderr, "%s refused: %s\n", nm.c_str(), why.c_str()); return 2; } rec = &X1[i]->rec; Qs[i] = X1[i]->Q; }
        else { bytes[i] = slurp(dir + "/" + nm + ".rv2"); X2[i].reset(new V2()); if (open_v2(bytes[i].data(), bytes[i].size(), C, *X2[i], why)) { fprintf(stderr, "%s refused: %s\n", nm.c_str(), why.c_str()); return 2; } rec = &X2[i]->rec; Qs[i] = X2[i]->Q; }
        for (auto& r : *rec) boff[i][r.hdr.substr(0, r.hdr.find_first_of(" \t"))] = r.boff; }
      std::vector<uint8_t>().swap(view.R); }
    const uint8_t* G = C.S.G; const uint64_t Rn = C.S.t[0].len;
    auto blk = [&](size_t i, uint64_t b, uint8_t* out, std::vector<RrOp>& ops, std::vector<uint8_t>& lit) -> int {
        if (X1[i]) { const V1& V = *X1[i]; const uint32_t blen = (uint32_t)std::min<uint64_t>(V.Q, V.nbases - b * V.Q);
            const int k = r3_decode_block(V.P + V.off[b], (uint32_t)(V.off[b + 1] - V.off[b]), V.T, G, Rn, blen, V.st[b], ops.data(), (uint32_t)ops.size(), lit.data(), (uint32_t)lit.size());
            if (k < 0) return 1; exec_fast(ops.data(), k, lit.data(), G, out); return 0; }
        return block_v2(*X2[i], C, b, out, ops, lit); };
    printf("OPEN\t%s\treference + %zu anchors + %zu archives\t%.3f s\n", label, an.size(), names.size(), now_s_() - t0); fflush(stdout);
    struct Scr { std::vector<RrOp> ops; std::vector<uint8_t> lit, blk; }; std::vector<Scr> sc(T); for (auto& s : sc) { s.ops.resize(RR_MAXOPS); s.lit.resize(16384); s.blk.resize(16384); }
    std::atomic<int> err{0}; int bad = 0;
    if (mode == "req") {
        std::vector<Req> q = read_reqs(argv[7]); std::vector<std::string> out(q.size()); uint64_t nb = 0; for (auto& r : q) nb += r.len;
        std::vector<size_t> qa(q.size()); std::vector<uint64_t> qs(q.size());
        for (size_t i = 0; i < q.size(); i++) { qa[i] = idx.at(q[i].sample); qs[i] = boff[qa[i]].at(q[i].ctg) + q[i].start; out[i].resize(q[i].len); }
        bad = protocol(label, WU, R, q.size(), nb, [&] { par_for(q.size(), T, [&](int t, size_t i) { const uint32_t Q = Qs[qa[i]]; Scr& s = sc[t]; const uint64_t a = qs[i], e = a + q[i].len;
                for (uint64_t b = a / Q; b <= (e - 1) / Q; b++) { if (blk(qa[i], b, s.blk.data(), s.ops, s.lit)) { err++; return; }
                    const uint64_t bs = b * Q, x0 = std::max(a, bs), x1 = std::min<uint64_t>(e, bs + Q); memcpy(&out[i][x0 - a], s.blk.data() + (x0 - bs), x1 - x0); } }); },
            [&] { size_t ok = 0; for (size_t i = 0; i < q.size(); i++) ok += !err && sha_upper(out[i].data(), q[i].len) == tr.req[q[i].id]; return ok; });
    } else {                                                                                  // each: one sample in memory at a time
        std::vector<std::string> S = split_c(argv[7]); double tot = 0; uint64_t nb = 0; size_t ok = 0;
        for (auto& smp : S) { const size_t i = idx.at(smp); const uint64_t nbs = X1[i] ? X1[i]->nbases : X2[i]->nbases, nbl = X1[i] ? X1[i]->nb : X2[i]->nb; const uint32_t Q = Qs[i];
            std::vector<uint8_t> A(nbs + 64); const double a0 = now_s_();
            par_for(nbl, T, [&](int t, size_t b) { if (blk(i, b, &A[b * Q], sc[t].ops, sc[t].lit)) err++; });
            const double dt = now_s_() - a0; tot += dt; nb += tr.bases[smp]; const bool good = !err && sha_upper((const char*)A.data(), tr.bases[smp]) == tr.sample_sha[smp]; ok += good;
            printf("EACH\t%s\t%s\t%s\t%.6f\t%llu bytes\t%s\n", label, smp.c_str(), X1[i] ? "anchor-v1" : "rv2", dt, (unsigned long long)tr.bases[smp], good ? "sha256 ok" : "sha256 DIFFERS"); fflush(stdout); }
        printf("RUN\t%s\t0\ttimed\t%.6f\t%zu items\t%llu bytes\tsha256 %zu/%zu\n", label, tot, S.size(), (unsigned long long)nb, ok, S.size()); bad = ok != S.size();
    }
    return bad ? 4 : 0;
}
// ---------------------------------------------------------------- byte ledger (H0 style) of v1 and rv2 archives
// ledger <t2t.fa> <v1 dir> <anchors|-> <rv2 dir> <names.txt> <mode v2|v1> <threads> - every archive parsed with an exact
// copy of its decoder's symbol order (r3_decode_stream for v1, r4_decode_stream for rv2) that charges each symbol
// log2(4096 / freq) bits and each raw bit 1 bit to its field; per archive: header (+ anchor table), meta (zstd, with its
// raw composition), block hashes, payload = sum of field bits / 8 + coder overhead (the residual: rANS state flush,
// renormalisation rounding, table quantisation is inside the field bits). mode v2: anchors = their v1 archives without
// the block-hash section (as counted), others rv2; mode v1: every name = its v1 q16k archive without block hashes.
enum { LF_LL, LF_LIT, LF_LITRAW, LF_KIND, LF_DELTA, LF_REP, LF_ABSDIR, LF_ABSTGT, LF_ABSPOS, LF_SELF, LF_FLIP, LF_LIFTTGT, LF_LIFTD,
       LF_LEN0, LF_LEN1, LF_LEN2, LF_LEN3, LF_LEN4, LF_LEN5, LF_LEN6, LF_N };
static const char* LF_NAME[LF_N] = {"literal_count", "literals", "literal_escape_raw", "kind", "delta_param", "rep_param", "abs_strand", "abs_target_id",
    "abs_position_raw32", "self_distance", "flip_param", "lift_target_id", "lift_delta", "len_CONT", "len_DELTA", "len_REP", "len_ABS", "len_SELF", "len_FLIP", "len_LIFT"};
struct AbsRec { uint8_t t; uint32_t pos; uint64_t dmin; };   // ABS study: target, local position, distance to the nearest cached locus in that target (~0 = none)
struct Led { double bits[LF_N] = {0}; uint64_t ev[7] = {0}; uint64_t blocks = 0, bad = 0; std::vector<AbsRec>* abs = nullptr;
    std::vector<uint32_t>* absp = nullptr; uint64_t basekind[7] = {0}; uint64_t litbases = 0; };
template <class Tab> static inline uint32_t led_sym(R3Dec* d, const Tab* T, int c, double* acc) {
    const uint32_t slot = d->x & (R3_M - 1); const uint32_t s = T->sym[c][slot]; const uint32_t f = T->freq[c][s]; if (!f) { d->bad = 1; return 0; }
    *acc += std::log2((double)R3_M / f); d->x = f * (d->x >> R3_PB) + slot - T->cum[c][s];
    while (d->x < R3_L) { if (d->p >= d->e) { d->bad = 1; return s; } d->x = (d->x << 8) | *d->p++; } return s; }
template <class Tab> static inline uint64_t led_val(R3Dec* d, const Tab* T, int c, double* acc) { const uint32_t s = led_sym(d, T, c, acc); const uint32_t nb = r3_nbits(s); *acc += nb; return r3_unbucket(s, r3_draw(d, nb)); }
// one block: V2 = rv2 contexts (rv2.h) with the target space; else v1 contexts (refrel3.h) against T2T only
template <bool V2, class Tab> static bool led_block(const uint8_t* src, uint32_t n, const Tab* T, const Space& S, uint64_t refn, uint32_t blen, R3Diag start, Led& L) {
    const int C_LL = 0, C_LITR = 4, C_LIT2 = 14, C_KIND = 31, C_DELTA = 47, C_REPK = 49, C_REPD = 50, C_DIR = 51, C_SELF = 52, C_LEN = 53;
    const int C_FLIPK = V2 ? 60 : 59, C_FLIPD = V2 ? 61 : 60;
    const uint8_t* ref = S.G; R3Dec d; r3_dinit(&d, src, n); if (d.bad) return false;
    R3Diag cache[4]; int nc = 1; cache[0] = start; int prevk = -1; uint32_t o = 0;
    auto kcl = [&](int k) { return V2 ? r4_kclass(k) : r3_kclass(k); };
    while (o < blen) {
        const uint64_t ll = led_val(&d, T, C_LL + kcl(prevk), &L.bits[LF_LL]); if (d.bad || ll > blen - o) return false;
        if (ll) { uint32_t p1 = 16, p2 = 16;
            for (uint32_t j = 0; j < ll; j++) { int c; if (ll <= 8) c = C_LITR + r3_refbase(ref, refn, cache[0], o + j) * 2 + (j == 0 ? 1 : 0); else c = C_LIT2 + (p1 == 16 ? 16 : (int)(p2 == 16 ? p1 : p2 * 4 + p1) % 16);
                const uint32_t s = led_sym(&d, T, c, &L.bits[LF_LIT]); if (s >= 5) { r3_draw(&d, 8); L.bits[LF_LITRAW] += 8; } p2 = p1; p1 = s < 4 ? s : 0; }
            o += (uint32_t)ll; L.litbases += ll; }
        if (o == blen) break;
        const int kind = (int)led_sym(&d, T, C_KIND + r3_llc(ll) * 4 + kcl(prevk), &L.bits[LF_KIND]);
        uint64_t p = 0; uint32_t dir = 0; int ki = 0; int64_t dd = 0; uint64_t dist = 0; int tt = 0;
        if (kind == K_CONT || kind == K_DELTA) { if (kind == K_DELTA) { const uint64_t z = led_val(&d, T, C_DELTA + (ll == 0 ? 0 : 1), &L.bits[LF_DELTA]); dd = (int64_t)(z >> 1) ^ -(int64_t)(z & 1); } }
        else if (kind == K_REP) { ki = 1 + (int)led_sym(&d, T, C_REPK, &L.bits[LF_REP]); const uint64_t z = led_val(&d, T, C_REPD, &L.bits[LF_REP]); dd = (int64_t)(z >> 1) ^ -(int64_t)(z & 1); if (ki >= nc) return false; }
        else if (kind == K_ABS) { dir = led_sym(&d, T, C_DIR, &L.bits[LF_ABSDIR]); if (V2) { tt = (int)led_sym(&d, T, R4_ABST, &L.bits[LF_ABSTGT]); if (tt >= (int)S.t.size()) return false; }
            p = r3_draw(&d, 32); L.bits[LF_ABSPOS] += 32;
            if (L.abs) { uint64_t dm = ~0ull; for (int k = 0; k < nc; k++) { const int64_t lo = r3_locus(cache[k], o); if (rv2_target(S, lo) == tt) { const uint64_t dd2 = (uint64_t)std::llabs(lo - (int64_t)(p + S.t[tt].base)); dm = std::min(dm, dd2); } }
                L.abs->push_back({(uint8_t)tt, (uint32_t)p, dm}); }
            if (L.absp) L.absp->push_back((uint32_t)p);
            p += S.t[tt].base; }
        else if (kind == K_SELF) dist = led_val(&d, T, C_SELF, &L.bits[LF_SELF]);
        else if (kind == K_FLIP) { ki = (int)led_sym(&d, T, C_FLIPK, &L.bits[LF_FLIP]); const uint64_t z = led_val(&d, T, C_FLIPD, &L.bits[LF_FLIP]); dd = (int64_t)(z >> 1) ^ -(int64_t)(z & 1); if (ki >= nc) return false; }
        else if (V2 && kind == K_LIFT) { tt = (int)led_sym(&d, T, R4_LIFTT, &L.bits[LF_LIFTTGT]); const uint64_t z = led_val(&d, T, R4_LIFTD + (ll == 0 ? 0 : 1), &L.bits[LF_LIFTD]); dd = (int64_t)(z >> 1) ^ -(int64_t)(z & 1); if (tt >= (int)S.t.size()) return false; }
        else return false;
        const uint64_t Ln = led_val(&d, T, C_LEN + kind, &L.bits[LF_LEN0 + kind]) + RR_MINL; if (d.bad || Ln > blen - o) return false;
        L.ev[kind]++; L.basekind[kind] += Ln;
        if (kind != K_SELF) {
            if (kind == K_FLIP) { const R3Diag g = cache[ki]; dir = g.dir ^ 1; const int64_t e = r3_locus(g, o) + dd; if (e < 0) return false; p = (uint64_t)e; }
            else if (kind == K_LIFT) { R3Diag ng; if (!rv2_lift(S, cache[0], o, tt, &ng)) return false; dir = ng.dir; const int64_t e = r3_exp(ng, o, Ln) + dd; if (e < 0) return false; p = (uint64_t)e; }
            else if (kind != K_ABS) { const R3Diag g = cache[kind == K_REP ? ki : 0]; dir = g.dir; const int64_t e = r3_exp(g, o, Ln) + dd; if (e < 0) return false; p = (uint64_t)e; }
            r3_push(cache, &nc, r3_diag(dir, p, o, Ln)); }
        else if (dist == 0 || dist > o) return false;
        o += (uint32_t)Ln; prevk = kind;
    }
    return d.x == R3_L && d.p == d.e;
}
static int rv2_cmd_ledger(char** argv) {
    const std::string v1dir = argv[3], dir = argv[5]; const std::vector<std::string> an = split_c(argv[4]); const std::vector<std::string> names = read_names(argv[6]);
    const bool v2 = !strcmp(argv[7], "v2"); const int T = atoi(argv[8]);
    Cohort C; load_cohort(C, argv[2], v1dir, v2 ? an : std::vector<std::string>(), T);
    Ref view; memcpy(view.sha, C.ref.sha, 32); view.name = C.ref.name; view.R.assign(C.G.begin(), C.G.begin() + C.S.t[0].len);
    // T2T-only space for the v1 archives (anchor offsets beyond T2T are never reached by v1 streams)
    Space S1; S1.G = C.S.G; S1.Gn = C.S.t[0].len; S1.t.resize(1); S1.t[0] = Tgt(); S1.t[0].base = 0; S1.t[0].len = C.S.t[0].len;
    struct Row { std::string name, role; uint64_t file = 0, hash = 0, counted = 0, header = 0, meta = 0, meta_raw = 0, payload = 0, nb = 0; uint64_t raw_ctg = 0, raw_case = 0, raw_model = 0, raw_blocks = 0; Led L; };
    std::vector<Row> rows(names.size()); std::atomic<size_t> next{0}; std::atomic<int> err{0};
    auto meta_parts = [](const uint8_t* mz, uint64_t mzn, uint64_t mr, int nctx, const int* alpha, uint64_t nb, Row& R) {
        std::vector<uint8_t> M(mr); ZSTD_decompress(M.data(), mr, mz, mzn); size_t i = 0;
        auto g32 = [&]() { uint32_t v; memcpy(&v, &M[i], 4); i += 4; return v; }; auto gl = [&]() { uint64_t v = 0; int sh = 0; for (;;) { uint8_t c = M[i++]; v |= (uint64_t)(c & 0x7F) << sh; if (!(c & 0x80)) break; sh += 7; } return v; };
        const uint32_t nr = g32(); for (uint32_t r = 0; r < nr; r++) { const uint32_t hl = g32(); i += hl + 12; } R.raw_ctg = i;
        uint64_t nl; memcpy(&nl, &M[i], 8); i += 8; for (uint64_t k = 0; k < nl; k++) { gl(); gl(); } R.raw_case = i - R.raw_ctg; const size_t a = i;
        for (int c = 0; c < nctx; c++) for (int s = 0; s < alpha[c]; s++) gl(); R.raw_model = i - a; R.raw_blocks = M.size() - i; (void)nb; };
    std::vector<std::thread> th; for (int t = 0; t < T; t++) th.emplace_back([&] { for (;;) { const size_t k = next.fetch_add(1); if (k >= names.size()) break;
        Row& R = rows[k]; R.name = names[k]; const bool isa = !v2 || std::find(an.begin(), an.end(), names[k]) != an.end();
        if (isa) { const std::vector<uint8_t> f = slurp(v1path(v1dir, names[k])); V1 X; std::string why; if (open_v1(f.data(), f.size(), view, X, why)) { err++; continue; }
            R.role = v2 ? "anchor(v1)" : "v1"; R.file = f.size(); R.hash = rd64(&f[112]); R.counted = R.file - R.hash; R.meta = rd64(&f[96]); R.meta_raw = rd64(&f[104]); R.payload = rd64(&f[120]); R.nb = X.nb;
            R.header = R.counted - R.meta - R.payload;
            meta_parts(f.data() + V1_HDR + 2 + (f[136] | f[137] << 8), R.meta, R.meta_raw, R3_NCTX, R3_ALPHA, X.nb, R);
            for (uint64_t b = 0; b < X.nb; b++) { const uint32_t blen = (uint32_t)std::min<uint64_t>(X.Q, X.nbases - b * X.Q);
                if (!led_block<false>(X.P + X.off[b], (uint32_t)(X.off[b + 1] - X.off[b]), X.T, S1, S1.Gn, blen, X.st[b], R.L)) R.L.bad++; R.L.blocks++; } }
        else { const std::vector<uint8_t> f = slurp(dir + "/" + names[k] + ".rv2"); V2 X; std::string why; if (open_v2(f.data(), f.size(), C, X, why)) { err++; continue; }
            R.role = "rv2"; R.file = R.counted = f.size(); R.meta = rd64(&f[96]); R.meta_raw = rd64(&f[104]); R.payload = rd64(&f[120]); R.nb = X.nb; R.header = R.counted - R.meta - R.payload;
            meta_parts(X.P - R.meta, R.meta, R.meta_raw, R4_NCTX, R4_ALPHA, X.nb, R);
            for (uint64_t b = 0; b < X.nb; b++) { const uint32_t blen = (uint32_t)std::min<uint64_t>(X.Q, X.nbases - b * X.Q);
                if (!led_block<true>(X.P + X.off[b], (uint32_t)(X.off[b + 1] - X.off[b]), X.T, C.S, C.S.Gn, blen, X.st[b], R.L)) R.L.bad++; R.L.blocks++; } } } });
    for (auto& x : th) x.join(); if (err) { fprintf(stderr, "%d archives refused\n", (int)err); return 2; }
    printf("LEDGER\tname\trole\tcounted\theader\tmeta\tpayload\tcoder_overhead\tblocks\tbad_blocks\tmeta_raw\traw_contigs\traw_case\traw_models\traw_blocktable");
    for (int f = 0; f < LF_N; f++) printf("\t%s", LF_NAME[f]); printf("\tev_CONT\tev_DELTA\tev_REP\tev_ABS\tev_SELF\tev_FLIP\tev_LIFT\n");
    Row tot; tot.name = "TOTAL"; tot.role = argv[7]; uint64_t bad = 0;
    for (auto& R : rows) { double fb = 0; for (int f = 0; f < LF_N; f++) fb += R.L.bits[f] / 8;
        printf("LEDGER\t%s\t%s\t%llu\t%llu\t%llu\t%llu\t%.0f\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu", R.name.c_str(), R.role.c_str(), (unsigned long long)R.counted, (unsigned long long)R.header, (unsigned long long)R.meta,
               (unsigned long long)R.payload, (double)R.payload - fb, (unsigned long long)R.L.blocks, (unsigned long long)R.L.bad, (unsigned long long)R.meta_raw, (unsigned long long)R.raw_ctg, (unsigned long long)R.raw_case, (unsigned long long)R.raw_model, (unsigned long long)R.raw_blocks);
        for (int f = 0; f < LF_N; f++) printf("\t%.0f", R.L.bits[f] / 8); for (int k = 0; k < 7; k++) printf("\t%llu", (unsigned long long)R.L.ev[k]); printf("\n");
        tot.counted += R.counted; tot.header += R.header; tot.meta += R.meta; tot.payload += R.payload; tot.L.blocks += R.L.blocks; bad += R.L.bad; tot.meta_raw += R.meta_raw;
        tot.raw_ctg += R.raw_ctg; tot.raw_case += R.raw_case; tot.raw_model += R.raw_model; tot.raw_blocks += R.raw_blocks; for (int f = 0; f < LF_N; f++) tot.L.bits[f] += R.L.bits[f]; for (int k = 0; k < 7; k++) tot.L.ev[k] += R.L.ev[k]; }
    double fb = 0; for (int f = 0; f < LF_N; f++) fb += tot.L.bits[f] / 8;
    printf("LEDGER\tTOTAL\t%s\t%llu\t%llu\t%llu\t%llu\t%.0f\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu", argv[7], (unsigned long long)tot.counted, (unsigned long long)tot.header, (unsigned long long)tot.meta, (unsigned long long)tot.payload,
           (double)tot.payload - fb, (unsigned long long)tot.L.blocks, (unsigned long long)bad, (unsigned long long)tot.meta_raw, (unsigned long long)tot.raw_ctg, (unsigned long long)tot.raw_case, (unsigned long long)tot.raw_model, (unsigned long long)tot.raw_blocks);
    for (int f = 0; f < LF_N; f++) printf("\t%.0f", tot.L.bits[f] / 8); for (int k = 0; k < 7; k++) printf("\t%llu", (unsigned long long)tot.L.ev[k]); printf("\n");
    return bad ? 3 : 0;
}
// ---------------------------------------------------------------- window law (model B of the hw-apex window-law job, protocol v1.1)
// lawtruth <t2t.fa> <v1 dir> <names.txt> <requests.tsv> > truth: every v1 archive decoded (all v1 checks, FASTA XXH3 ==
//   header), then SAMPLE name bases XXH3(upper bases) and REQ W id SHA-256(upper bases of the request)
// law <v1|v1h|rv2> <t2t.fa> <v1 dir> <anchors|-> <rv2 dir> <names.txt> <requests.tsv> <truth> - one thread, in process:
//   D_Q = canonical bases of the cohort / median of 9 full sequential decodes (3 warm-ups; output buffer allocated
//   outside the timer; XXH3 of every decoded assembly checked after its interval); windows: per request one timed fetch
//   (contig lookup by name + the touched blocks decoded + copy), SHA-256 checked after the timer; p50 nearest rank;
//   c0 = p50(W = 1) - Q / D_Q; predicted t(W) = c0 + (W + Q - 1) / D_Q; error = 100 (p50 - pred) / pred; PASS if
//   c0 >= 0 and |error| <= 20 % at W = 1024, 8192, 65536. v1 = v1 archives, block XXH3 not checked; v1h = checked;
//   rv2 = rv2 archives of the non-anchors + v1 archives of the anchors (block XXH3 not checked).
struct LawReq { int W; std::string id, asmb, ctg; uint64_t s, e; };
static std::vector<LawReq> law_reqs(const char* p) { std::vector<LawReq> v; for (auto& r : read_tsv(p)) if (r.size() >= 6 && r[0] != "W") v.push_back({atoi(r[0].c_str()), r[1], r[2], r[3], strtoull(r[4].c_str(), 0, 10), strtoull(r[5].c_str(), 0, 10)}); return v; }
static int rv2_cmd_lawtruth(char** argv) {
    Ref ref = load_ref(argv[2], 1, false); const std::string v1dir = argv[3]; const std::vector<std::string> names = read_names(argv[4]); const std::vector<LawReq> q = law_reqs(argv[5]);
    for (auto& nm : names) { const std::vector<uint8_t> f = slurp(v1path(v1dir, nm)); V1 X; std::string why, fa;
        if (open_v1(f.data(), f.size(), ref, X, why) || full_v1(X, ref, fa, true, 16)) { fprintf(stderr, "%s: v1 decode failed\n", nm.c_str()); return 2; }
        std::string().swap(fa); std::vector<uint8_t> A(X.nbases + 64); std::vector<RrOp> ops(RR_MAXOPS); std::vector<uint8_t> lit(X.Q);
        for (uint64_t b = 0; b < X.nb; b++) if (block_v1(X, ref, b, &A[b * X.Q], ops, lit, true)) { fprintf(stderr, "%s block %llu\n", nm.c_str(), (unsigned long long)b); return 3; }
        printf("SAMPLE\t%s\t%llu\t%016llx\n", nm.c_str(), (unsigned long long)X.nbases, (unsigned long long)XXH3_64bits(A.data(), X.nbases));
        std::map<std::string, uint64_t> bo; for (auto& r : X.rec) bo[r.hdr.substr(0, r.hdr.find_first_of(" \t"))] = r.boff;
        for (auto& r : q) if (r.asmb == nm) { unsigned char md[32]; const uint64_t o = bo.at(r.ctg); SHA256(&A[o + r.s], r.e - r.s, md); printf("REQ\t%d\t%s\t%s\n", r.W, r.id.c_str(), sha_hex(md).c_str()); }
        fflush(stdout); }
    return 0;
}
// wintruth <t2t.fa> <v1 dir> <requests.tsv (id sample contig start length)> <names map: manifest <tab> set name> > truth:
// REQ id SHA-256 of the upper-case window, from the v1 archives with every block XXH3 checked (block_v1 verify) -
// truth for diagnostic request sets without a full decode of the cohort
static int rv2_cmd_wintruth(char** argv) {
    Ref ref = load_ref(argv[2], 1, false); const std::string v1dir = argv[3]; const std::vector<Req> q = read_reqs(argv[4]);
    std::map<std::string, std::string> man; for (auto& v : read_tsv(argv[5])) man[v[1]] = v[0];
    std::map<std::string, std::unique_ptr<V1>> X; std::map<std::string, std::vector<uint8_t>> by; std::vector<RrOp> ops(RR_MAXOPS); std::vector<uint8_t> lit(16384), blk(16384);
    for (auto& r : q) { const std::string nm = man.count(r.sample) ? man[r.sample] : r.sample;
        if (!X.count(nm)) { by[nm] = slurp(v1path(v1dir, nm)); X[nm].reset(new V1()); std::string why; if (open_v1(by[nm].data(), by[nm].size(), ref, *X[nm], why)) { fprintf(stderr, "%s refused\n", nm.c_str()); return 2; } }
        const V1& V = *X[nm]; uint64_t bo = ~0ull; for (auto& c : V.rec) if (c.hdr.substr(0, c.hdr.find_first_of(" \t")) == r.ctg) bo = c.boff;
        if (bo == ~0ull) { fprintf(stderr, "no contig %s\n", r.ctg.c_str()); return 2; }
        const uint64_t a = bo + r.start, e = a + r.len; std::string out;
        for (uint64_t b = a / V.Q; b <= (e - 1) / V.Q; b++) { if (block_v1(V, ref, b, blk.data(), ops, lit, true)) { fprintf(stderr, "block failed\n"); return 3; }
            const uint64_t bs = b * V.Q, x0 = std::max(a, bs), x1 = std::min<uint64_t>(e, bs + V.Q); out.append((const char*)blk.data() + (x0 - bs), x1 - x0); }
        unsigned char md[32]; SHA256((const unsigned char*)out.data(), out.size(), md); printf("REQ\t%s\t%s\n", r.id.c_str(), sha_hex(md).c_str()); }
    return 0;
}
// absstudy <t2t.fa> <v1 dir> <anchors> <rv2 dir> <names.txt> - lever 2 offline: every ABS event of the rv2 archives
// (target, 32-bit local position, distance to the nearest cached locus of the same target); empirical code lengths of
// candidate models (static per cohort, an upper estimate of a per-archive table): raw 32 bits; top 8 bits by target +
// 24 raw; top 12 bits by target + 20 raw; bucketed distance to the nearest same-target cached locus (+ sign) when one
// exists, else top-12 model, with a 1-bit choice
static int rv2_cmd_absstudy(char** argv) {
    const std::vector<std::string> an = split_c(argv[4]); Cohort C; load_cohort(C, argv[2], argv[3], an, 16);
    std::vector<AbsRec> rec; uint64_t nfiles = 0;
    for (auto& nm : read_names(argv[6])) { if (std::find(an.begin(), an.end(), nm) != an.end()) continue;
        const std::vector<uint8_t> f = slurp(std::string(argv[5]) + "/" + nm + ".rv2"); V2 X; std::string why; if (open_v2(f.data(), f.size(), C, X, why)) return 2; nfiles++;
        Led L; L.abs = &rec; for (uint64_t b = 0; b < X.nb; b++) { const uint32_t blen = (uint32_t)std::min<uint64_t>(X.Q, X.nbases - b * X.Q);
            if (!led_block<true>(X.P + X.off[b], (uint32_t)(X.off[b + 1] - X.off[b]), X.T, C.S, C.S.Gn, blen, X.st[b], L)) return 3; } }
    const double n = (double)rec.size(); auto H = [&](std::map<uint64_t, uint64_t>& m) { double h = 0; for (auto& kv : m) h -= kv.second * std::log2(kv.second / n); return h; };
    std::map<uint64_t, uint64_t> m8, m12, mb; double rawd = 0; uint64_t withc = 0;
    for (auto& r : rec) { m8[(uint64_t)r.t << 8 | r.pos >> 24]++; m12[(uint64_t)r.t << 12 | r.pos >> 20]++; }
    double hyb = 0; std::map<uint64_t, uint64_t> mh;   // hybrid: choice bit + (bucket of distance + sign) or top-12
    for (auto& r : rec) { if (r.dmin != ~0ull && r.dmin < (1ull << 26)) { uint32_t nb2; uint64_t ex; const uint32_t sy = r3_bucket(r.dmin, &nb2, &ex); mh[1ull << 40 | sy]++; hyb += nb2 + 1; withc++; }
        else { mh[(uint64_t)r.t << 12 | r.pos >> 20]++; hyb += 20; } }
    const double h8 = H(m8) + 24 * n, h12 = H(m12) + 20 * n, hh = H(mh) + hyb;
    printf("ABSSTUDY\tarchives %llu\tABS events %.0f\traw32 %.0f B\ttop8|t+24 %.0f B\ttop12|t+20 %.0f B\thybrid(near cached locus <2^26: %llu) %.0f B\n", (unsigned long long)nfiles, n, n * 4, h8 / 8, h12 / 8, (unsigned long long)withc, hh / 8);
    (void)rawd; return 0;
}
static int rv2_cmd_law(char** argv) {
    const std::string var = argv[2], v1dir = argv[4], dir = argv[6]; const std::vector<std::string> an = var == "rv2" ? split_c(argv[5]) : std::vector<std::string>();
    const std::vector<std::string> names = read_names(argv[7]); const std::vector<LawReq> q = law_reqs(argv[8]); const bool hv = var == "v1h";
    std::map<std::string, std::pair<uint64_t, uint64_t>> ts; std::map<std::pair<int, std::string>, std::string> tr;
    for (auto& r : read_tsv(argv[9])) { if (r[0] == "SAMPLE") ts[r[1]] = {strtoull(r[2].c_str(), 0, 10), strtoull(r[3].c_str(), 0, 16)}; else if (r[0] == "REQ") tr[{atoi(r[1].c_str()), r[2]}] = r[3]; }
    const double t0 = now_s(); Cohort C; load_cohort(C, argv[3], v1dir, an, 16);
    Ref view; memcpy(view.sha, C.ref.sha, 32); view.name = C.ref.name; view.R.assign(C.G.begin(), C.G.begin() + C.S.t[0].len);
    const size_t N = names.size(); std::vector<std::vector<uint8_t>> bytes(N); std::vector<std::unique_ptr<V1>> X1(N); std::vector<std::unique_ptr<V2>> X2(N);
    std::vector<std::vector<std::string>> rn(N); std::vector<const std::vector<Rec>*> rec(N); std::map<std::string, size_t> idx; uint64_t total = 0, maxb = 0;
    for (size_t i = 0; i < N; i++) { const std::string& nm = names[i]; idx[nm] = i; std::string why;
        if (var != "rv2" || std::find(an.begin(), an.end(), nm) != an.end()) { bytes[i] = slurp(v1path(v1dir, nm)); X1[i].reset(new V1()); if (open_v1(bytes[i].data(), bytes[i].size(), view, *X1[i], why)) { fprintf(stderr, "%s refused: %s\n", nm.c_str(), why.c_str()); return 2; } rec[i] = &X1[i]->rec; }
        else { bytes[i] = slurp(dir + "/" + nm + ".rv2"); X2[i].reset(new V2()); if (open_v2(bytes[i].data(), bytes[i].size(), C, *X2[i], why)) { fprintf(stderr, "%s refused: %s\n", nm.c_str(), why.c_str()); return 2; } rec[i] = &X2[i]->rec; }
        for (auto& r : *rec[i]) rn[i].push_back(r.hdr.substr(0, r.hdr.find_first_of(" \t")));
        const uint64_t nb = X1[i] ? X1[i]->nbases : X2[i]->nbases; total += nb; maxb = std::max(maxb, nb);
        if (!ts.count(nm) || ts[nm].first != nb) { fprintf(stderr, "%s: no truth / size differs\n", nm.c_str()); return 2; } }
    std::vector<uint8_t>().swap(view.R);
    const uint8_t* G = C.S.G; const uint64_t Rn = C.S.t[0].len; const uint32_t Q = 16384;
    std::vector<RrOp> ops(RR_MAXOPS); std::vector<uint8_t> lit(Q), blk(Q);
    auto block = [&](size_t i, uint64_t b, uint8_t* out) -> int {
        if (X1[i]) { const V1& V = *X1[i];
            const uint32_t blen = (uint32_t)std::min<uint64_t>(V.Q, V.nbases - b * V.Q);
            const int k = r3_decode_block(V.P + V.off[b], (uint32_t)(V.off[b + 1] - V.off[b]), V.T, G, Rn, blen, V.st[b], ops.data(), (uint32_t)ops.size(), lit.data(), (uint32_t)lit.size());
            if (k < 0) return 1; exec_fast(ops.data(), k, lit.data(), G, out); if (hv && V.hashes && XXH3_64bits(out, blen) != rd64(V.hashes + 8 * b)) return 2; return 0; }
        return block_v2(*X2[i], C, b, out, ops, lit); };
    for (size_t i = 0; i < N; i++) if ((X1[i] ? X1[i]->Q : X2[i]->Q) != Q) { fprintf(stderr, "Q differs\n"); return 2; }
    printf("LAW\tvariant %s\tanchors %zu\tassemblies %zu\tcanonical bases %llu\topen %.1f s\n", var.c_str(), an.size(), N, (unsigned long long)total, now_s() - t0); fflush(stdout);
    auto ns = [] { return (uint64_t)std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now().time_since_epoch()).count(); };
    // D_Q
    std::vector<uint8_t> A(maxb + 64); std::vector<uint64_t> dq;
    for (int it = 0; it < 12; it++) { uint64_t sum = 0; bool ok = true;
        for (size_t i = 0; i < N; i++) { const uint64_t nb = X1[i] ? X1[i]->nb : X2[i]->nb, nbs = X1[i] ? X1[i]->nbases : X2[i]->nbases;
            const uint64_t a = ns(); for (uint64_t b = 0; b < nb; b++) if (block(i, b, &A[b * Q])) { ok = false; break; } const uint64_t z = ns(); sum += z - a;
            if (!ok || XXH3_64bits(A.data(), nbs) != ts[names[i]].second) { printf("DQ\t%d\tFAILED %s\n", it, names[i].c_str()); return 4; } }
        printf("DQ\t%d\t%s\t%llu ns\tall %zu XXH3 ok\n", it, it < 3 ? "warmup" : "measured", (unsigned long long)sum, N); fflush(stdout); if (it >= 3) dq.push_back(sum); }
    std::vector<uint64_t> s9 = dq; std::sort(s9.begin(), s9.end()); const double med = (double)s9[4], DQ = total / (med * 1e-9);
    printf("DQ_RESULT\tmedian %.0f ns\tmin %llu\tmax %llu\tD_Q %.6f Gbases/s\n", med, (unsigned long long)s9[0], (unsigned long long)s9[8], DQ / 1e9);
    // windows
    std::vector<uint8_t> dst(65536 + 64); double c0 = 0; int fails = 0; std::string raw;
    for (int W : {1, 1024, 8192, 65536}) { std::vector<uint64_t> lat; int bad = 0;
        for (auto& r : q) { if (r.W != W) continue; const size_t i = idx.at(r.asmb);
            const uint64_t a = ns();
            size_t ri = 0; while (ri < rn[i].size() && rn[i][ri] != r.ctg) ri++;
            if (ri == rn[i].size() || r.e > (*rec[i])[ri].len) { bad++; continue; }
            const uint64_t lo = (*rec[i])[ri].boff + r.s, hi = (*rec[i])[ri].boff + r.e; size_t n = 0; bool ok = true;
            for (uint64_t b = lo / Q; b <= (hi - 1) / Q; b++) { if (block(i, b, blk.data())) { ok = false; break; } const uint64_t bs = b * Q, x0 = std::max(lo, bs), x1 = std::min<uint64_t>(hi, bs + Q); memcpy(&dst[n], blk.data() + x0 - bs, x1 - x0); n += x1 - x0; }
            const uint64_t z = ns();
            unsigned char md[32]; SHA256(dst.data(), n, md); if (!ok || n != r.e - r.s || sha_hex(md) != tr[{W, r.id}]) { bad++; continue; }
            lat.push_back(z - a); raw += std::to_string(W) + "\t" + r.id + "\t" + std::to_string(z - a) + "\n"; }
        std::vector<uint64_t> s = lat; std::sort(s.begin(), s.end()); auto nr = [&](int p) { return s.empty() ? 0.0 : s[(p * s.size() + 99) / 100 - 1] / 1000.0; };
        const double p50 = nr(50); if (W == 1) c0 = p50 - Q / DQ * 1e6;
        const double pred = c0 + (W + Q - 1) / DQ * 1e6, err = 100 * (p50 - pred) / pred; const bool verdict = W != 1;
        const bool pass = bad == 0 && (!verdict || (c0 >= 0 && std::fabs(err) <= 20.0)); if (!pass) fails++;
        printf("WIN\tW %d\tverified %zu\tbad %d\tp50 %.3f us\tp95 %.3f\tp99 %.3f\tpred %.3f us\terror %+.2f %%\t%s\n", W, lat.size(), bad, p50, nr(95), nr(99), pred, err,
               W == 1 ? "c0-probe" : pass ? "PASS" : "FAIL"); fflush(stdout); }
    printf("LAW_RESULT\tvariant %s\tQ %u\tD_Q %.6f Gbases/s\tc0 %.3f us\tverdict %s\n", var.c_str(), Q, DQ / 1e9, c0, fails ? "FAIL" : "PASS");
    if (getenv("LAW_RAW")) spit(getenv("LAW_RAW"), raw.data(), raw.size());
    return fails ? 5 : 0;
}
static int rv2_cmd_absv1(char** argv) {
    Ref ref = load_ref(argv[2], 1, false); const std::string v1dir = argv[3];
    const std::vector<std::string> names = read_names(argv[4]); const int T = atoi(argv[5]); const std::string outdir = argv[6];
    Space S1; S1.G = ref.R.data(); S1.Gn = ref.R.size(); S1.t.resize(1); S1.t[0] = Tgt(); S1.t[0].base = 0; S1.t[0].len = ref.R.size();
    printf("name\tbases\tblocks\tbad\tev_CONT\tev_DELTA\tev_REP\tev_ABS\tev_SELF\tev_FLIP\tb_CONT\tb_DELTA\tb_REP\tb_ABS\tb_SELF\tb_FLIP\tb_literal\tabs_bits\tfile_bits\n");
    for (auto& nm : names) {
        const std::vector<uint8_t> f = slurp(v1path(v1dir, nm)); V1 X; std::string why;
        if (open_v1(f.data(), f.size(), ref, X, why)) { fprintf(stderr, "%s refused: %s\n", nm.c_str(), why.c_str()); return 2; }
        std::vector<std::vector<uint32_t>> pth(T); std::vector<Led> Lth(T);
        std::atomic<uint64_t> next{0}; std::vector<std::thread> th;
        for (int t = 0; t < T; t++) th.emplace_back([&, t] { Lth[t].absp = &pth[t];
            for (;;) { const uint64_t b = next.fetch_add(1); if (b >= X.nb) break;
                const uint32_t blen = (uint32_t)std::min<uint64_t>(X.Q, X.nbases - b * X.Q);
                if (!led_block<false>(X.P + X.off[b], (uint32_t)(X.off[b + 1] - X.off[b]), X.T, S1, S1.Gn, blen, X.st[b], Lth[t])) Lth[t].bad++;
                Lth[t].blocks++; } });
        for (auto& x : th) x.join();
        Led A; std::vector<uint32_t> all;
        for (auto& L : Lth) { A.blocks += L.blocks; A.bad += L.bad; A.litbases += L.litbases;
            for (int k = 0; k < 7; k++) { A.ev[k] += L.ev[k]; A.basekind[k] += L.basekind[k]; }
            for (int i = 0; i < LF_N; i++) A.bits[i] += L.bits[i];
            all.insert(all.end(), pth[&L - Lth.data()].begin(), pth[&L - Lth.data()].end()); }
        const std::string op = outdir + "/" + nm + ".abs32";
        FILE* o = fopen(op.c_str(), "wb"); if (!o) { fprintf(stderr, "cannot write %s\n", op.c_str()); return 3; }
        fwrite(all.data(), 4, all.size(), o); fclose(o);
        printf("%s\t%llu\t%llu\t%llu", nm.c_str(), (unsigned long long)X.nbases, (unsigned long long)A.blocks, (unsigned long long)A.bad);
        for (int k = 0; k < 6; k++) printf("\t%llu", (unsigned long long)A.ev[k]);
        for (int k = 0; k < 6; k++) printf("\t%llu", (unsigned long long)A.basekind[k]);
        printf("\t%llu\t%.0f\t%llu\n", (unsigned long long)A.litbases, A.bits[LF_ABSPOS], (unsigned long long)(f.size() - rd64(&f[112])) * 8);
        fflush(stdout);
    }
    return 0;
}
// locate <t2t.fa> <v1 dir> <names.txt> <threads> <n positions> <seed> <verify n> - (assembly, position) ->
// (T2T position, strand) for a base copied from the reference, SELF -> the position its chain resolves to, LIT -> none.
// Uses only the op list of the block the position falls in: kind 1 reference, 3 reference reverse complement,
// 2 self copy (source is earlier in the same block), 0 literal.
struct LocRes { uint64_t ref = 0; int cls = 0; uint32_t hops = 0; uint8_t strand = 0; };   // cls 1 ref, 2 self-resolved, 0 literal, -1 error
static inline int loc_in_block(const RrOp* ops, int nk, uint32_t o, LocRes* r) {
    uint32_t hops = 0; uint32_t cur = o; uint32_t flip = 0;
    for (;;) {
        int found = -1;
        for (int i = 0; i < nk; i++) if (cur >= ops[i].dst && cur < ops[i].dst + ops[i].len) { found = i; break; }
        if (found < 0) { r->cls = -1; return 1; }
        const RrOp& e = ops[found]; const uint32_t d = cur - e.dst;
        if (e.kind == 0) { r->cls = 0; r->hops = hops; return 0; }
        if (e.kind == 1) { r->ref = e.src + d; r->strand = (uint8_t)flip; r->cls = hops ? 2 : 1; r->hops = hops; return 0; }
        if (e.kind == 3) { r->ref = e.src + e.len - 1 - d; r->strand = (uint8_t)(flip ^ 1); r->cls = hops ? 2 : 1; r->hops = hops; return 0; }
        cur = (uint32_t)e.src + d; hops++;                                  // self copy: source is block-local
        if (hops > 100000) { r->cls = -1; return 1; }
    }
}
static int rv2_cmd_locate(char** argv) {
    Ref ref = load_ref(argv[2], 1, false); const std::string v1dir = argv[3];
    const std::vector<std::string> names = read_names(argv[4]);
    const int T = atoi(argv[5]); const uint64_t NP = strtoull(argv[6], 0, 10);
    const uint64_t seed = strtoull(argv[7], 0, 10); const uint64_t VN = strtoull(argv[8], 0, 10);
    const size_t NA = names.size();
    std::vector<std::vector<uint8_t>> files(NA); std::vector<std::unique_ptr<V1>> X(NA); std::vector<uint64_t> cum(NA + 1, 0);
    for (size_t i = 0; i < NA; i++) { files[i] = slurp(v1path(v1dir, names[i])); X[i].reset(new V1()); std::string why;
        if (open_v1(files[i].data(), files[i].size(), ref, *X[i], why)) { fprintf(stderr, "%s refused: %s\n", names[i].c_str(), why.c_str()); return 2; }
        cum[i + 1] = cum[i] + X[i]->nbases; }
    const uint64_t TOT = cum[NA];
    printf("LOCATE\tassemblies %zu\tbases %llu\tpositions %llu\tthreads %d\n", NA, (unsigned long long)TOT, (unsigned long long)NP, T);
    std::vector<uint64_t> pos(NP); std::vector<uint32_t> asm_(NP);
    { std::mt19937_64 g(seed); for (uint64_t i = 0; i < NP; i++) { const uint64_t x = g() % TOT;
        const size_t a = (size_t)(std::upper_bound(cum.begin(), cum.end(), x) - cum.begin() - 1); asm_[i] = (uint32_t)a; pos[i] = x - cum[a]; } }
    std::vector<uint64_t> cls(4, 0), hopsum(1, 0); std::vector<uint64_t> clsth(4 * T, 0); std::vector<uint64_t> hopth(T, 0), errth(T, 0);
    std::atomic<uint64_t> next{0}; const double t0 = now_s();
    std::vector<std::thread> th;
    for (int t = 0; t < T; t++) th.emplace_back([&, t] { std::vector<RrOp> ops(RR_MAXOPS); std::vector<uint8_t> lit(16384);
        for (;;) { const uint64_t i = next.fetch_add(256); if (i >= NP) break;
            for (uint64_t k = i; k < std::min(i + 256, NP); k++) {
                const V1& V = *X[asm_[k]]; const uint64_t b = pos[k] / V.Q; const uint32_t o = (uint32_t)(pos[k] % V.Q);
                const uint32_t blen = (uint32_t)std::min<uint64_t>(V.Q, V.nbases - b * V.Q);
                const int nk = r3_decode_block(V.P + V.off[b], (uint32_t)(V.off[b + 1] - V.off[b]), V.T, ref.R.data(), ref.R.size(), blen, V.st[b], ops.data(), (uint32_t)ops.size(), lit.data(), (uint32_t)lit.size());
                if (nk < 0) { errth[t]++; continue; }
                LocRes r; if (loc_in_block(ops.data(), nk, o, &r)) { errth[t]++; continue; }
                clsth[4 * t + (r.cls < 0 ? 3 : r.cls)]++; hopth[t] += r.hops; } } });
    for (auto& x : th) x.join();
    const double dt = now_s() - t0;
    uint64_t err = 0, hops = 0; for (int t = 0; t < T; t++) { for (int c = 0; c < 4; c++) cls[c] += clsth[4 * t + c]; err += errth[t]; hops += hopth[t]; }
    printf("RATE\tthreads %d\t%.3f s\t%.0f positions/s\n", T, dt, NP / dt);
    printf("CLASS\tliteral %llu\tref_direct %llu\tself_resolved %llu\tunresolved %llu\terrors %llu\tself hops total %llu\n",
           (unsigned long long)cls[0], (unsigned long long)cls[1], (unsigned long long)cls[2], (unsigned long long)cls[3], (unsigned long long)err, (unsigned long long)hops);
    // verification: the assembly base at the position must equal the reference base at the located coordinate
    uint64_t ok = 0, bad = 0, skipped = 0; const uint64_t VM = std::min(VN, NP);
    { std::vector<RrOp> ops(RR_MAXOPS); std::vector<uint8_t> lit(16384); std::vector<uint8_t> blk(16384);
      for (uint64_t k = 0; k < VM; k++) {
        const V1& V = *X[asm_[k]]; const uint64_t b = pos[k] / V.Q; const uint32_t o = (uint32_t)(pos[k] % V.Q);
        const uint32_t blen = (uint32_t)std::min<uint64_t>(V.Q, V.nbases - b * V.Q);
        const int nk = r3_decode_block(V.P + V.off[b], (uint32_t)(V.off[b + 1] - V.off[b]), V.T, ref.R.data(), ref.R.size(), blen, V.st[b], ops.data(), (uint32_t)ops.size(), lit.data(), (uint32_t)lit.size());
        if (nk < 0) { bad++; continue; }
        exec_fast(ops.data(), nk, lit.data(), ref.R.data(), blk.data());
        LocRes r; if (loc_in_block(ops.data(), nk, o, &r) || r.cls <= 0) { skipped++; continue; }
        const uint8_t asmb = blk[o]; const uint8_t refb = r.strand ? rr_comp(ref.R[r.ref]) : ref.R[r.ref];
        if (asmb == refb) ok++; else { if (bad < 5) printf("MISMATCH\t%s\tpos %llu\tT2T %llu\tstrand %c\tassembly %c\treference %c\n",
            names[asm_[k]].c_str(), (unsigned long long)pos[k], (unsigned long long)r.ref, r.strand ? '-' : '+', asmb, refb); bad++; } } }
    printf("VERIFY\tchecked %llu\tidentical %llu\tdiffer %llu\tskipped (literal) %llu\n",
           (unsigned long long)(ok + bad), (unsigned long long)ok, (unsigned long long)bad, (unsigned long long)skipped);
    return bad ? 4 : 0;
}
// tilemap <t2t.fa> <v1 dir> <names.txt> <W> <threads> <out.tsv> - tile every contig of every assembly into
// non-overlapping windows of W bases and give each window its interval in the reference: locate the first and the
// last base of the window, and report the hull [lo, hi) plus the strand and whether both ends are reference-derived.
// Windows never cross a contig. A window whose two ends land far apart (a rearrangement or a long insertion) is
// reported as it is, with its span, so the caller can decide what to do with it.
static int rv2_cmd_tilemap(char** argv) {
    Ref ref = load_ref(argv[2], 1, false); const std::string v1dir = argv[3];
    const std::vector<std::string> names = read_names(argv[4]);
    const uint32_t W = (uint32_t)strtoul(argv[5], 0, 10); const int T = atoi(argv[6]);
    FILE* o = fopen(argv[7], "w"); if (!o) { fprintf(stderr, "cannot write %s\n", argv[7]); return 2; }
    fprintf(o, "assembly\tcontig\tstart\tref_lo\tref_hi\tstrand\tcls_start\tcls_end\n");
    uint64_t nwin = 0, nboth = 0, nlit = 0, nspan = 0; const double t0 = now_s();
    for (size_t ai = 0; ai < names.size(); ai++) {
        const std::vector<uint8_t> f = slurp(v1path(v1dir, names[ai])); V1 X; std::string why;
        if (open_v1(f.data(), f.size(), ref, X, why)) { fprintf(stderr, "%s refused: %s\n", names[ai].c_str(), why.c_str()); return 3; }
        struct Row { uint32_t ci; uint64_t st; uint64_t lo, hi; char strand; int cs, ce; };
        std::vector<std::vector<Row>> out(T);
        std::vector<std::pair<uint32_t, uint64_t>> jobs;           // contig index, window start inside the contig
        for (uint32_t c = 0; c < X.rec.size(); c++) for (uint64_t s = 0; s + W <= X.rec[c].len; s += W) jobs.push_back({c, s});
        std::atomic<size_t> next{0}; std::vector<std::thread> th;
        for (int t = 0; t < T; t++) th.emplace_back([&, t] { std::vector<RrOp> ops(RR_MAXOPS); std::vector<uint8_t> lit(16384);
            auto one = [&](uint64_t pos, LocRes* r) {
                const uint64_t b = pos / X.Q; const uint32_t off = (uint32_t)(pos % X.Q);
                const uint32_t blen = (uint32_t)std::min<uint64_t>(X.Q, X.nbases - b * X.Q);
                const int nk = r3_decode_block(X.P + X.off[b], (uint32_t)(X.off[b + 1] - X.off[b]), X.T, ref.R.data(), ref.R.size(), blen, X.st[b], ops.data(), (uint32_t)ops.size(), lit.data(), (uint32_t)lit.size());
                if (nk < 0) { r->cls = -1; return; }
                if (loc_in_block(ops.data(), nk, off, r)) r->cls = -1; };
            for (;;) { const size_t j = next.fetch_add(64); if (j >= jobs.size()) break;
                for (size_t k = j; k < std::min(j + 64, jobs.size()); k++) {
                    const uint64_t base = X.rec[jobs[k].first].boff + jobs[k].second;
                    LocRes a, b2; one(base, &a); one(base + W - 1, &b2);
                    Row r; r.ci = jobs[k].first; r.st = jobs[k].second; r.cs = a.cls; r.ce = b2.cls;
                    r.strand = (a.cls > 0 && a.strand) ? '-' : '+';
                    if (a.cls > 0 && b2.cls > 0) { r.lo = std::min(a.ref, b2.ref); r.hi = std::max(a.ref, b2.ref) + 1; }
                    else { r.lo = 0; r.hi = 0; }
                    out[t].push_back(r); } } });
        for (auto& x : th) x.join();
        for (int t = 0; t < T; t++) for (const Row& r : out[t]) {
            nwin++; if (r.hi > r.lo) { nboth++; const uint64_t sp = r.hi - r.lo; if (sp > (uint64_t)W + W / 10 || sp < (uint64_t)W - W / 10) nspan++; } else nlit++;
            fprintf(o, "%s\t%s\t%llu\t%llu\t%llu\t%c\t%d\t%d\n", names[ai].c_str(), X.rec[r.ci].hdr.substr(0, X.rec[r.ci].hdr.find_first_of(" \t")).c_str(),
                    (unsigned long long)r.st, (unsigned long long)r.lo, (unsigned long long)r.hi, r.strand, r.cs, r.ce); }
        fprintf(stderr, "%s: %zu windows\n", names[ai].c_str(), jobs.size());
    }
    fclose(o);
    printf("TILEMAP\tW %u\twindows %llu\tboth ends located %llu\tan end in a literal %llu\tspan off by >10%% %llu\t%.1f s\t%.0f windows/s\n",
           W, (unsigned long long)nwin, (unsigned long long)nboth, (unsigned long long)nlit, (unsigned long long)nspan, now_s() - t0, nwin / (now_s() - t0));
    return 0;
}
// leak <t2t.fa> <v1 dir> <names.txt> <W> <threads> <people|locus> <seed> <chunk Mb> - leakage by locus.
// Every contig is tiled into non-overlapping windows of W bases; each window gets its interval in the reference
// (locate the first and the last base, take the hull) and is used only if the hull is within 10 % of W, so a window
// that crosses a rearrangement is not pretended to be one locus.
//   people: train / test split by INDIVIDUAL (the sample name without the haplotype suffix), random, about half each.
//   locus:  the reference is cut into chunks of <chunk Mb>, each chunk goes to train or test at random; a window
//           belongs to the side of the chunk its interval starts in.
// Then: the share of test windows whose interval is covered by the union of the train windows' intervals by >= 50 %
// and >= 90 %. Two passes over the cohort; the train side is kept as a merged interval list, which stays small.
struct Iv { uint64_t lo, hi; };
static void iv_merge_into(std::vector<Iv>& acc, std::vector<Iv>& add) {
    if (add.empty()) return;
    std::sort(add.begin(), add.end(), [](const Iv& a, const Iv& b) { return a.lo < b.lo; });
    std::vector<Iv> out; out.reserve(acc.size() + add.size());
    size_t i = 0, j = 0;
    auto push = [&](const Iv& v) { if (!out.empty() && v.lo <= out.back().hi) { if (v.hi > out.back().hi) out.back().hi = v.hi; } else out.push_back(v); };
    while (i < acc.size() || j < add.size()) {
        if (j >= add.size() || (i < acc.size() && acc[i].lo <= add[j].lo)) push(acc[i++]); else push(add[j++]);
    }
    acc.swap(out);
}
static std::string individual_of(const std::string& nm) {        // y1_HG00438.1 -> HG00438 ; r2_HG00097_hap1_... -> HG00097
    size_t a = nm.find('_'); std::string t = a == std::string::npos ? nm : nm.substr(a + 1);
    size_t b = t.find_first_of("._"); return b == std::string::npos ? t : t.substr(0, b);
}
static int rv2_cmd_leak(char** argv) {
    Ref ref = load_ref(argv[2], 1, false); const std::string v1dir = argv[3];
    const std::vector<std::string> names = read_names(argv[4]);
    const uint32_t W = (uint32_t)strtoul(argv[5], 0, 10); const int T = atoi(argv[6]);
    const std::string mode = argv[7]; const uint64_t seed = strtoull(argv[8], 0, 10);
    const uint64_t CH = strtoull(argv[9], 0, 10) * 1000000ull;
    std::mt19937_64 g(seed);
    std::vector<char> side(names.size(), 0);                      // 1 = train (people mode)
    std::vector<uint64_t> train_chunks;                           // locus mode: sorted chunk ids in train
    if (mode == "people") {
        std::vector<std::string> inds; for (auto& n : names) { const std::string i2 = individual_of(n); if (std::find(inds.begin(), inds.end(), i2) == inds.end()) inds.push_back(i2); }
        std::vector<size_t> ord(inds.size()); for (size_t i = 0; i < ord.size(); i++) ord[i] = i;
        for (size_t i = ord.size(); i > 1; i--) std::swap(ord[i - 1], ord[g() % i]);
        std::vector<char> tr(inds.size(), 0); for (size_t i = 0; i < inds.size() / 2; i++) tr[ord[i]] = 1;
        for (size_t k = 0; k < names.size(); k++) side[k] = tr[(size_t)(std::find(inds.begin(), inds.end(), individual_of(names[k])) - inds.begin())];
        size_t ntr = 0; for (char c : side) ntr += c;
        printf("LEAK\tmode people\tindividuals %zu\ttrain individuals %zu\ttrain assemblies %zu of %zu\n", inds.size(), inds.size() / 2, ntr, names.size());
    } else {
        const uint64_t nch = ref.R.size() / CH + 1;
        for (uint64_t c = 0; c < nch; c++) if (g() & 1) train_chunks.push_back(c);
        printf("LEAK\tmode locus\tchunk %llu Mb\tchunks %llu\ttrain chunks %zu\n", (unsigned long long)(CH / 1000000), (unsigned long long)nch, train_chunks.size());
    }
    auto in_train_chunk = [&](uint64_t lo) { const uint64_t c = lo / CH; return std::binary_search(train_chunks.begin(), train_chunks.end(), c); };

    std::vector<Iv> train;                                        // merged
    uint64_t nwin = 0, nuse = 0, nskip_span = 0, nskip_lit = 0;
    uint64_t ntest = 0, ge50 = 0, ge90 = 0, zero = 0; long double covsum = 0;
    const double t0 = now_s();
    for (int pass = 0; pass < 2; pass++) {
        if (pass == 1 && mode == "people") { /* both passes needed in locus mode only, but people mode is cheap too */ }
        for (size_t ai = 0; ai < names.size(); ai++) {
            if (mode == "people" && ((pass == 0) != (side[ai] == 1))) continue;    // pass 0 = train assemblies, pass 1 = test
            const std::vector<uint8_t> f = slurp(v1path(v1dir, names[ai])); V1 X; std::string why;
            if (open_v1(f.data(), f.size(), ref, X, why)) { fprintf(stderr, "%s refused: %s\n", names[ai].c_str(), why.c_str()); return 3; }
            std::vector<std::pair<uint32_t, uint64_t>> jobs;
            for (uint32_t c = 0; c < X.rec.size(); c++) for (uint64_t s = 0; s + W <= X.rec[c].len; s += W) jobs.push_back({c, s});
            std::vector<std::vector<Iv>> mine(T);
            std::vector<uint64_t> cnt(6 * T, 0);                  // per thread: win, use, skip_span, skip_lit, test, (unused)
            std::vector<uint64_t> acc(4 * T, 0);                  // per thread: test, ge50, ge90, zero
            std::vector<long double> cs(T, 0);
            std::atomic<size_t> next{0}; std::vector<std::thread> th;
            for (int t = 0; t < T; t++) th.emplace_back([&, t] { std::vector<RrOp> ops(RR_MAXOPS); std::vector<uint8_t> lit(16384);
                auto one = [&](uint64_t pos, LocRes* r) {
                    const uint64_t b = pos / X.Q; const uint32_t off = (uint32_t)(pos % X.Q);
                    const uint32_t blen = (uint32_t)std::min<uint64_t>(X.Q, X.nbases - b * X.Q);
                    const int nk = r3_decode_block(X.P + X.off[b], (uint32_t)(X.off[b + 1] - X.off[b]), X.T, ref.R.data(), ref.R.size(), blen, X.st[b], ops.data(), (uint32_t)ops.size(), lit.data(), (uint32_t)lit.size());
                    if (nk < 0) { r->cls = -1; return; }
                    if (loc_in_block(ops.data(), nk, off, r)) r->cls = -1; };
                for (;;) { const size_t j = next.fetch_add(64); if (j >= jobs.size()) break;
                    for (size_t k = j; k < std::min(j + 64, jobs.size()); k++) {
                        const uint64_t base = X.rec[jobs[k].first].boff + jobs[k].second;
                        LocRes a, b2; one(base, &a); one(base + W - 1, &b2);
                        cnt[6 * t + 0]++;
                        if (a.cls <= 0 || b2.cls <= 0) { cnt[6 * t + 3]++; continue; }
                        const uint64_t lo = std::min(a.ref, b2.ref), hi = std::max(a.ref, b2.ref) + 1;
                        const uint64_t span = hi - lo;
                        if (span > (uint64_t)W + W / 10 || span < (uint64_t)W - W / 10) { cnt[6 * t + 2]++; continue; }
                        cnt[6 * t + 1]++;
                        const bool is_train = mode == "people" ? side[ai] == 1 : in_train_chunk(lo);
                        if (pass == 0) { if (is_train) mine[t].push_back({lo, hi}); continue; }
                        if (is_train) continue;                   // pass 1 counts the test side only
                        uint64_t cov = 0;                         // covered by the merged train intervals
                        size_t i2 = (size_t)(std::upper_bound(train.begin(), train.end(), lo, [](uint64_t v, const Iv& x) { return v < x.hi; }) - train.begin());
                        for (; i2 < train.size() && train[i2].lo < hi; i2++) {
                            const uint64_t x = std::max(train[i2].lo, lo), y = std::min(train[i2].hi, hi);
                            if (x < y) cov += y - x; }
                        const double frac = (double)cov / (double)span;
                        acc[4 * t + 0]++; if (frac >= 0.5) acc[4 * t + 1]++; if (frac >= 0.9) acc[4 * t + 2]++; if (cov == 0) acc[4 * t + 3]++;
                        cs[t] += frac; } } });
            for (auto& x : th) x.join();
            for (int t = 0; t < T; t++) { nwin += cnt[6 * t + 0]; nuse += cnt[6 * t + 1]; nskip_span += cnt[6 * t + 2]; nskip_lit += cnt[6 * t + 3];
                ntest += acc[4 * t + 0]; ge50 += acc[4 * t + 1]; ge90 += acc[4 * t + 2]; zero += acc[4 * t + 3]; covsum += cs[t]; }
            if (pass == 0) for (int t = 0; t < T; t++) iv_merge_into(train, mine[t]);
        }
        if (pass == 0) { uint64_t tb = 0; for (const Iv& v : train) tb += v.hi - v.lo;
            printf("LEAK\ttrain merged intervals %zu\treference bases covered %llu (%.1f %% of T2T)\n", train.size(), (unsigned long long)tb, 100.0 * tb / ref.R.size()); fflush(stdout); }
    }
    const double dt = now_s() - t0;
    printf("LEAK\tW %u\twindows seen %llu\tusable %llu\tskipped: span %llu, literal end %llu\n",
           W, (unsigned long long)nwin, (unsigned long long)nuse, (unsigned long long)nskip_span, (unsigned long long)nskip_lit);
    printf("LEAK\ttest windows %llu\t>=50%% overlap %llu (%.4f %%)\t>=90%% overlap %llu (%.4f %%)\tzero overlap %llu (%.4f %%)\tmean overlap %.4f %%\n",
           (unsigned long long)ntest, (unsigned long long)ge50, 100.0 * ge50 / std::max<uint64_t>(1, ntest),
           (unsigned long long)ge90, 100.0 * ge90 / std::max<uint64_t>(1, ntest),
           (unsigned long long)zero, 100.0 * zero / std::max<uint64_t>(1, ntest), 100.0 * (double)(covsum / std::max<uint64_t>(1, ntest)));
    printf("LEAK\twall %.1f s\twindows/s %.0f\tthreads %d\n", dt, nwin / dt, T);
    return 0;
}
int main(int argc, char** argv) {
    for (int c = 0; c < 256; c++) COMP[c] = rr_comp((uint8_t)c);
    if (argc >= 7 && !strcmp(argv[1], "absstudy")) return rv2_cmd_absstudy(argv);
    if (argc >= 6 && !strcmp(argv[1], "wintruth")) return rv2_cmd_wintruth(argv);
    if (argc >= 10 && !strcmp(argv[1], "leak")) return rv2_cmd_leak(argv);
    if (argc >= 8 && !strcmp(argv[1], "tilemap")) return rv2_cmd_tilemap(argv);
    if (argc >= 9 && !strcmp(argv[1], "locate")) return rv2_cmd_locate(argv);
    if (argc >= 7 && !strcmp(argv[1], "absv1")) return rv2_cmd_absv1(argv);
    if (argc >= 9 && !strcmp(argv[1], "ledger")) return rv2_cmd_ledger(argv);
    if (argc >= 6 && !strcmp(argv[1], "lawtruth")) return rv2_cmd_lawtruth(argv);
    if (argc >= 10 && !strcmp(argv[1], "law")) return rv2_cmd_law(argv);
    if (argc >= 9 && !strcmp(argv[1], "judge")) return rv2_cmd_judge(argv);
    if (argc >= 3 && !strcmp(argv[1], "info")) return rv2_cmd_info(argv);
    if (argc >= 3 && !strcmp(argv[1], "bench")) return rv2_cmd_bench(argc, argv);
    fprintf(stderr, "usage: rv2 judge|bench|info|law|lawtruth|wintruth|ledger|absstudy ... (see the header)\n"); return 1;
}
