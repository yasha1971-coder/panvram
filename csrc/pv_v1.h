// SPDX-License-Identifier: MIT
// pv_v1.h - refrel3 v1 reader (FORMAT.md), host side: the reference (FASTA -> upper-case base stream + SHA-256), an
// archive opened with every check of FORMAT.md section 5 before any decode, the CPU block decoder, full decode -> FASTA.
// Ported from aceapex research/refrel @ 5b6d5ce: open_v1 / block_v1 / full_v1 (refrel3v1.cpp), build_tab (refrel3.cpp),
// exec_fast (refrel.cpp), parse_fasta (refrel_io.h). Changes: errors are exceptions (no exit), the archive's arrays can
// live in caller-owned pooled memory (Archive::attach), Q is the archive's.
#pragma once
#include "refrel3.h"
#include "pv_compact.h"
#include "rr_sha256.h"
#define XXH_INLINE_ALL
#include "xxhash.h"
#include <zstd.h>
#include <algorithm>
#include <cstdio>
#include <cstring>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace pv {

static const char V1_MAGIC[8] = {'R', 'F', 'R', 'L', '3', 'V', '1', 0};
static const uint32_t V1_VERSION = 1, V1_HDR = 136, F_BLOCKHASH = 1;
static inline bool q_ok(uint32_t q) { return q == 1024 || q == 2048 || q == 4096 || q == 16384; }
static inline uint64_t rd64(const uint8_t* p) { uint64_t v; memcpy(&v, p, 8); return v; }
static inline uint32_t rd32(const uint8_t* p) { uint32_t v; memcpy(&v, p, 4); return v; }
static inline std::string hex(const uint8_t* p, size_t n) { static const char* d = "0123456789abcdef"; std::string s; for (size_t i = 0; i < n; i++) { s += d[p[i] >> 4]; s += d[p[i] & 15]; } return s; }

struct Err : std::runtime_error { using std::runtime_error::runtime_error; };

static std::vector<uint8_t> slurp(const std::string& p) {
    FILE* f = fopen(p.c_str(), "rb"); if (!f) throw Err(p + ": cannot open");
    fseek(f, 0, SEEK_END); const long s = ftell(f); fseek(f, 0, SEEK_SET); std::vector<uint8_t> v((size_t)s);
    const size_t r = fread(v.data(), 1, v.size(), f); fclose(f); if (r != v.size()) throw Err(p + ": short read");
    return v;
}
// FASTA -> its bases in place: header lines dropped, line ends dropped, a..z upper-cased (the decoded reference of
// FORMAT.md 1). Returns the base count; buf[0, count) holds them.
static uint64_t fasta_bases_inplace(uint8_t* buf, uint64_t n, const std::string& path) {
    if (!n || buf[0] != '>') throw Err(path + ": not FASTA (no header at byte 0)");
    uint64_t i = 0, w = 0;
    while (i < n) {
        if (buf[i] == '>') { const uint8_t* e = (const uint8_t*)memchr(buf + i, '\n', n - i); if (!e) throw Err(path + ": header without line end"); i = (uint64_t)(e - buf) + 1; continue; }
        const uint8_t* e = (const uint8_t*)memchr(buf + i, '\n', n - i); const uint64_t le = e ? (uint64_t)(e - buf) : n;
        if (!e) throw Err(path + ": no final line end");
        for (uint64_t k = i; k < le; k++) { const uint8_t c = buf[k]; buf[w++] = (c >= 'a' && c <= 'z') ? (uint8_t)(c - 32) : c; }
        i = le + 1;
    }
    return w;
}

// ---------------------------------------------------------------- tables, execution
static inline void build_tab(R3Tab& T) {
    for (int c = 0; c < R3_NCTX; c++) { uint32_t a = 0; for (int s = 0; s < R3_AB; s++) { T.cum[c][s] = (uint16_t)a; for (uint32_t k = 0; k < T.freq[c][s]; k++) T.sym[c][a + k] = (uint8_t)s; a += T.freq[c][s]; } T.cum[c][R3_AB] = (uint16_t)a; }
}
struct CompTable { uint8_t t[256]; CompTable() { for (int c = 0; c < 256; c++) t[c] = rr_comp((uint8_t)c); } };
static const CompTable COMP;
static inline void ref_copy(uint8_t* o, const uint8_t* ref, uint64_t src, uint32_t len) { memcpy(o, ref + src, len); }
static inline void ref_copy(uint8_t* o, const PackedRef& ref, uint64_t src, uint32_t len) { for (uint32_t j = 0; j < len; j++) o[j] = ref[src + j]; }
static inline uint8_t ref_at(const uint8_t* ref, uint64_t q) { return ref[q]; }
static inline uint8_t ref_at(const PackedRef& ref, uint64_t q) { return ref[q]; }
template <class RP> static inline void exec_ops(const RrOp* ops, int n, const uint8_t* lit, const RP& ref, uint8_t* out) {
    for (int k = 0; k < n; k++) { const RrOp& q = ops[k];
        if (q.kind == 0) memcpy(out + q.dst, lit + q.src, q.len);
        else if (q.kind == 1) ref_copy(out + q.dst, ref, q.src, q.len);
        else if (q.kind == 2) { const uint32_t dist = q.dst - (uint32_t)q.src;
            if (dist >= q.len) memcpy(out + q.dst, out + q.src, q.len); else for (uint32_t j = 0; j < q.len; j++) out[q.dst + j] = out[q.src + j]; }
        else { uint8_t* o = out + q.dst; const uint64_t r = q.src + q.len - 1; for (uint32_t j = 0; j < q.len; j++) o[j] = COMP.t[ref_at(ref, r - j)]; } }
}
// output transforms shared by the CPU path and (as device functions) the kernels: complement keeping case, tokens
static inline uint8_t comp_case(uint8_t c) { switch (c) { case 'A': return 'T'; case 'C': return 'G'; case 'G': return 'C'; case 'T': return 'A';
    case 'a': return 't'; case 'c': return 'g'; case 'g': return 'c'; case 't': return 'a'; default: return c; } }
static inline uint8_t token(uint8_t c) { switch (c | 0x20) { case 'a': return 0; case 'c': return 1; case 'g': return 2; case 't': return 3; default: return 4; } }

// ---------------------------------------------------------------- archive
struct Contig { std::string hdr, name; uint64_t len, boff; uint32_t lw; };
struct Archive {
    std::string path, refname; uint32_t Q = 0, flags = 0; uint64_t nbases = 0, nb = 0, fasta_xxh = 0, bases_xxh = 0, file_bytes = 0, payload_bytes = 0;
    std::vector<Contig> rec; std::vector<uint64_t> low_s, low_l;                 // lower-case runs, sorted, disjoint
    std::vector<uint64_t> hashes;                                               // per block XXH3 (flag bit 0), else empty
    // owned while opening, then moved into pooled memory by the caller (attach)
    std::unique_ptr<R3Tab> T_own; std::vector<uint32_t> off_own; std::vector<int64_t> st_own; std::vector<uint8_t> P_own;
    const R3Tab* T = nullptr; const uint32_t* off = nullptr; const int64_t* st = nullptr; const uint8_t* P = nullptr;
    void attach(const R3Tab* t, const uint32_t* o, const int64_t* s, const uint8_t* p) {
        T = t; off = o; st = s; P = p; T_own.reset(); std::vector<uint32_t>().swap(off_own); std::vector<int64_t>().swap(st_own); std::vector<uint8_t>().swap(P_own); }
    R3Diag state(uint64_t b) const { R3Diag g; g.c = (uint64_t)(st[b] >> 1); g.dir = (uint32_t)(st[b] & 1); return g; }
    // compact block table (opt-in): offsets and states from CompactTab; off / st unused then
    const CompactTab* ct = nullptr; uint64_t gb0 = 0, g0 = 0;
    void locate(uint64_t b, uint32_t* o, uint32_t* l, R3Diag* g) const {
        if (ct) { ct->get(gb0, g0, b, o, l, g); return; }
        *o = off[b]; *l = off[b + 1] - off[b]; *g = state(b); }
};
static inline int64_t pack_state(R3Diag g) { return (int64_t)((uint64_t)g.c << 1) | (int64_t)(g.dir & 1); }

// open_v1 of refrel3v1.cpp: every check before any decode; throws Err("<path>: <why>")
static void open_v1(const std::vector<uint8_t>& file, const std::string& path, const uint8_t ref_sha[32], uint64_t ref_n, Archive& X) {
    const uint8_t* a = file.data(); const size_t n = file.size();
    auto fail = [&](const char* w) { throw Err(path + ": refused: " + w); };
    if (n < V1_HDR + 2) fail("short file");
    if (memcmp(a, V1_MAGIC, 8)) fail("magic (not refrel3 v1)");
    if (rd32(a + 8) != V1_VERSION) fail("version");
    X.Q = rd32(a + 12); if (!q_ok(X.Q)) fail("block size");
    X.flags = rd32(a + 16); if (X.flags & ~F_BLOCKHASH) fail("unknown flags"); if (rd32(a + 20)) fail("reserved");
    X.nbases = rd64(a + 64); X.nb = rd64(a + 72); X.fasta_xxh = rd64(a + 80); X.bases_xxh = rd64(a + 88);
    const uint64_t mz = rd64(a + 96), mr = rd64(a + 104), hs = rd64(a + 112), pl = rd64(a + 120);
    const uint16_t nl = (uint16_t)(a[136] | a[137] << 8);
    const uint64_t pstart = (uint64_t)V1_HDR + 2 + nl + mz + hs;
    if (mz > n || hs > n || pl > n || pstart > n || pstart + pl != n) fail("section sizes");
    if (X.nb != (X.nbases + X.Q - 1) / X.Q || (X.nbases && !X.nb)) fail("block count");
    if (hs != ((X.flags & F_BLOCKHASH) ? X.nb * 8 : 0)) fail("hash section size");
    { std::vector<uint8_t> h(a, a + pstart); memset(&h[128], 0, 8); if (XXH3_64bits(h.data(), h.size()) != rd64(a + 128)) fail("header XXH3"); }
    if (memcmp(a + 24, ref_sha, 32) || rd64(a + 56) != ref_n) fail("reference SHA-256 / size differs (wrong reference)");
    X.refname.assign((const char*)a + 138, nl);
    const uint8_t* mzp = a + V1_HDR + 2 + nl;
    if (mr > (1ull << 32) || ZSTD_getFrameContentSize(mzp, mz) != mr) fail("meta frame");
    if (ZSTD_findFrameCompressedSize(mzp, mz) != mz) fail("meta frame count (exactly one zstd frame)");
    std::vector<uint8_t> M(mr); if (ZSTD_decompress(M.data(), mr, mzp, mz) != mr) fail("meta decompress");
    size_t i = 0; bool bad = false;
    auto need = [&](size_t k) { if (i + k > M.size()) { bad = true; return false; } return true; };
    auto g64 = [&]() -> uint64_t { if (!need(8)) return 0; uint64_t v; memcpy(&v, &M[i], 8); i += 8; return v; };
    auto g32 = [&]() -> uint32_t { if (!need(4)) return 0; uint32_t v; memcpy(&v, &M[i], 4); i += 4; return v; };
    auto gl = [&]() -> uint64_t { uint64_t v = 0; int sh = 0; for (;;) { if (!need(1) || sh > 63) { bad = true; return 0; } uint8_t c = M[i++]; v |= (uint64_t)(c & 0x7F) << sh; if (!(c & 0x80)) break; sh += 7; } return v; };
    const uint32_t nr = g32(); uint64_t bo = 0; if (nr > M.size()) fail("records");
    for (uint32_t r = 0; r < nr && !bad; r++) { Contig q; const uint32_t hl = g32(); if (!need(hl)) break; q.hdr.assign((const char*)&M[i], hl); i += hl; q.len = g64(); q.lw = g32(); if (!q.lw && q.len) bad = true;
        q.name = q.hdr.substr(0, q.hdr.find_first_of(" \t")); q.boff = bo; bo += q.len; X.rec.push_back(q); }
    if (bad || bo != X.nbases) fail("contig table");
    const uint64_t nlw = g64(); uint64_t lst = 0; if (nlw > X.nbases) fail("case runs");
    for (uint64_t k = 0; k < nlw && !bad; k++) { const uint64_t g = gl(), l = gl(); if (lst + g + l > X.nbases) { bad = true; break; } X.low_s.push_back(lst + g); X.low_l.push_back(l); lst += g + l; }
    if (bad) fail("case runs");
    X.T_own.reset(new R3Tab()); memset(X.T_own.get(), 0, sizeof(R3Tab));
    for (int c = 0; c < R3_NCTX && !bad; c++) { uint64_t sum = 0; for (int s = 0; s < R3_ALPHA[c]; s++) { const uint64_t f = gl(); if (f > R3_M) bad = true; X.T_own->freq[c][s] = (uint16_t)f; sum += f; } if (sum && sum != R3_M) bad = true; }
    if (bad) fail("model tables");
    build_tab(*X.T_own);
    if (pl >= (1ull << 32)) fail("payload over 4 GiB (pooled offsets are 32-bit)");
    X.off_own.assign(X.nb + 1, 0); X.st_own.assign(X.nb + 1, 0); uint64_t prevc = 0; uint32_t prevd = 0, acc = 0;
    for (uint64_t b = 0; b < X.nb && !bad; b++) { const uint64_t len = gl(); if (len > pl - acc) { bad = true; break; } acc += (uint32_t)len; X.off_own[b + 1] = acc; const uint64_t v = gl();
        const int64_t pred = b == 0 ? 0 : (int64_t)(prevd ? prevc - X.Q : prevc + X.Q); const uint64_t q = v >> 1;
        R3Diag g; g.dir = (uint32_t)(v & 1); g.c = (uint64_t)(pred + ((int64_t)(q >> 1) ^ -(int64_t)(q & 1))); X.st_own[b] = pack_state(g); prevc = g.c; prevd = g.dir; }
    if (bad || i != M.size() || acc != pl) fail("block table");
    if (X.flags & F_BLOCKHASH) { X.hashes.resize(X.nb); memcpy(X.hashes.data(), a + V1_HDR + 2 + nl + mz, X.nb * 8); }
    X.P_own.assign(a + pstart, a + n); X.payload_bytes = pl; X.file_bytes = n; X.path = path;
    X.T = X.T_own.get(); X.off = X.off_own.data(); X.st = X.st_own.data(); X.P = X.P_own.data();
}

// block b (upper case) into out; 0 ok, 1 decode error, 2 block XXH3 mismatch (verify and hashes present)
struct Scratch { std::vector<RrOp> ops; std::vector<uint8_t> lit, blk; explicit Scratch(uint32_t Q) : ops(RR_MAXOPS), lit(Q), blk(Q) {} };
template <class RP> static inline int decode_block(const Archive& X, const RP& ref, uint64_t ref_n, uint64_t b, uint8_t* out, Scratch& S, bool verify) {
    const uint32_t blen = (uint32_t)std::min<uint64_t>(X.Q, X.nbases - b * X.Q);
    uint32_t o, l; R3Diag g; X.locate(b, &o, &l, &g);
    const int k = r3_decode_block(X.P + o, l, X.T, ref, ref_n, blen, g, S.ops.data(), (uint32_t)S.ops.size(), S.lit.data(), (uint32_t)S.lit.size());
    if (k < 0) return 1;
    exec_ops(S.ops.data(), k, S.lit.data(), ref, out);
    if (verify && !X.hashes.empty() && XXH3_64bits(out, blen) != X.hashes[b]) return 2;
    return 0;
}
// first lower-case run ending after x
static inline size_t first_run(const Archive& X, uint64_t x) {
    size_t lo = 0, hi = X.low_s.size(); while (lo < hi) { const size_t m = (lo + hi) / 2; if (X.low_s[m] + X.low_l[m] <= x) lo = m + 1; else hi = m; } return lo;
}
// bases [s, s + W) of the stream into out, case applied if asked; 0 ok, else decode_block's code
template <class RP> static inline int window_cpu(const Archive& X, const RP& ref, uint64_t ref_n, uint64_t s, uint64_t W, uint8_t* out, Scratch& S, bool apply_case, bool verify) {
    if (!W) return 0;
    for (uint64_t b = s / X.Q; b <= (s + W - 1) / X.Q; b++) {
        const int e = decode_block(X, ref, ref_n, b, S.blk.data(), S, verify); if (e) return e;
        const uint64_t bs = b * X.Q, x0 = std::max(s, bs), x1 = std::min<uint64_t>(s + W, bs + X.Q);
        memcpy(out + (x0 - s), S.blk.data() + (x0 - bs), x1 - x0);
    }
    if (apply_case) for (size_t r = first_run(X, s); r < X.low_s.size() && X.low_s[r] < s + W; r++) {
        const uint64_t x0 = std::max(s, X.low_s[r]), x1 = std::min(s + W, X.low_s[r] + X.low_l[r]); for (uint64_t x = x0; x < x1; x++) out[x - s] |= 0x20; }
    return 0;
}

}  // namespace pv
