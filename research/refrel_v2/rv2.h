// rv2.h - refrel v2 PROTOTYPE (H4 "cohort anchors"; NOT a frozen format): the v1 grammar of refrel3.h over a TARGET
// SPACE = T2T (target 0) followed by K anchor assemblies (targets 1..K, their decoded upper-case bases), plus:
//   - ABS carries a target id (symbol, context R4_ABST) before its 32-bit position, which is local to that target;
//   - a new kind LIFT (6): target id (symbol, R4_LIFTT) + zigzag delta (R4_LIFTD): the copy goes to the position of
//     target t that corresponds, through the anchors' own copy lists, to the base just before the current offset on the
//     current diagonal (cache[0]); see rv2_lift below;
//   - every copy must lie inside one target.
// Depth 1: anchors are refrel3 v1 archives (T2T only); their decoded bases and their reference-copy lists (the "lift
// tables") are resident in the decoder. With K = 0 the stream is the v1 grammar with 66 instead of 61 contexts.
#pragma once
#include "refrel3.h"
#include <vector>
#include <algorithm>

enum { R4_LL = 0, R4_LITR = 4, R4_LIT2 = 14, R4_KIND = 31, R4_DELTA = 47, R4_REPK = 49, R4_REPD = 50, R4_DIR = 51, R4_SELF = 52, R4_LEN = 53,
       R4_FLIPK = 60, R4_FLIPD = 61, R4_ABST = 62, R4_LIFTT = 63, R4_LIFTD = 64, R4_NCTX = 66 };
enum { K_LIFT = 6 };
#define R4_MAXT 16                                                     // targets incl. T2T (target symbol alphabet)
static const int R4_ALPHA[R4_NCTX] = {
    R3_AB, R3_AB, R3_AB, R3_AB,
    6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
    6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,                 // KIND: 7 kinds
    R3_AB, R3_AB, 3, R3_AB, 2, R3_AB,
    R3_AB, R3_AB, R3_AB, R3_AB, R3_AB, R3_AB, R3_AB,                 // LEN x 7 kinds
    4, R3_AB,
    R4_MAXT, R4_MAXT, R3_AB, R3_AB };                                // ABST, LIFTT, LIFTD (ll == 0, else)
static inline int r4_kclass(int k) { return k < 0 ? 0 : k == K_CONT ? 1 : (k == K_DELTA || k == K_REP || k == K_FLIP || k == K_LIFT) ? 2 : 3; }
struct R4Tab { uint16_t freq[R4_NCTX][R3_AB]; uint16_t cum[R4_NCTX][R3_AB + 1]; uint8_t sym[R4_NCTX][R3_M]; };
static inline uint32_t r4_dsym(R3Dec* d, const R4Tab* T, int c) {
    const uint32_t slot = d->x & (R3_M - 1); const uint32_t s = T->sym[c][slot];
    const uint32_t f = T->freq[c][s]; if (!f) { d->bad = 1; return 0; }
    d->x = f * (d->x >> R3_PB) + slot - T->cum[c][s];
    while (d->x < R3_L) { if (d->p >= d->e) { d->bad = 1; return s; } d->x = (d->x << 8) | *d->p++; }
    return s;
}
static inline uint64_t r4_dval(R3Dec* d, const R4Tab* T, int c) { const uint32_t s = r4_dsym(d, T, c); return r3_unbucket(s, r3_draw(d, r3_nbits(s))); }

// ---------------------------------------------------------------- target space and lift tables
// LiftE: one reference copy of an anchor (from its v1 stream, in stream order): anchor[a, a+L) = T2T[p, p+L) (dir 0)
// or the reverse complement of T2T[p, p+L) (dir 1). a is local to the anchor.
struct LiftE { uint64_t p, a; uint32_t L, dir; };
struct Tgt { uint64_t base = 0, len = 0; std::vector<LiftE> byA; std::vector<uint32_t> byP; };   // byP: indices sorted by (p, a)
struct Space { const uint8_t* G = nullptr; uint64_t Gn = 0; std::vector<Tgt> t; };
static inline int rv2_target(const Space& S, int64_t x) {
    if (x < 0) return -1; for (int k = (int)S.t.size() - 1; k >= 0; k--) if ((uint64_t)x >= S.t[k].base) return (uint64_t)x < S.t[k].base + S.t[k].len ? k : -1; return -1; }
static inline bool rv2_inside(const Space& S, uint64_t p, uint64_t L) { const int k = rv2_target(S, (int64_t)p); return k >= 0 && L <= S.t[k].base + S.t[k].len - p; }
// anchor-local y -> T2T x (flip ^= entry dir): the last entry with a <= y (first entry if none), extrapolated
static inline bool rv2_a2t(const Tgt& T, int64_t y, int64_t* x, int* flip) {
    const std::vector<LiftE>& v = T.byA; if (v.empty()) return false;
    size_t lo = 0, hi = v.size(); while (lo < hi) { const size_t m = (lo + hi) / 2; if ((int64_t)v[m].a <= y) lo = m + 1; else hi = m; }
    const LiftE& e = v[lo ? lo - 1 : 0];
    if (!e.dir) *x = (int64_t)e.p + (y - (int64_t)e.a); else { *x = (int64_t)e.p + e.L - 1 - (y - (int64_t)e.a); *flip ^= 1; }
    return true;
}
// T2T x -> anchor-local y: among the (up to 16) last entries in (p, a) order with p <= x, the longest one containing x
// (first in backward order on ties); if none contains x, the last entry with p <= x (first entry if none); extrapolated
static inline bool rv2_t2a(const Tgt& T, int64_t x, int64_t* y, int* flip) {
    const std::vector<uint32_t>& ix = T.byP; if (ix.empty()) return false;
    size_t lo = 0, hi = ix.size(); while (lo < hi) { const size_t m = (lo + hi) / 2; if ((int64_t)T.byA[ix[m]].p <= x) lo = m + 1; else hi = m; }
    size_t pick = lo ? lo - 1 : 0; uint32_t bestL = 0;
    for (size_t k = lo, n = 0; k > 0 && n < 16; k--, n++) { const LiftE& e = T.byA[ix[k - 1]]; if (x < (int64_t)(e.p + e.L) && e.L > bestL) { bestL = e.L; pick = k - 1; } }
    const LiftE& e = T.byA[ix[pick]];
    if (!e.dir) *y = (int64_t)e.a + (x - (int64_t)e.p); else { *y = (int64_t)e.a + ((int64_t)e.p + e.L - 1 - x); *flip ^= 1; }
    return true;
}
// LIFT prediction: the base before block offset o on diagonal g (global xg) mapped into target tt -> the new diagonal
// (its strand = g's strand ^ the flips of the mapping), aligned so that offset o - 1 sits on the mapped base.
static inline bool rv2_lift(const Space& S, R3Diag g, uint32_t o, int tt, R3Diag* ng) {
    const int64_t xg = g.dir ? (int64_t)g.c - (int64_t)o + 1 : (int64_t)g.c + (int64_t)o - 1;
    const int s = rv2_target(S, xg); if (s < 0 || tt < 0 || tt >= (int)S.t.size()) return false;
    int64_t yg; int flip = 0;
    if (s == tt) yg = xg;
    else { int64_t x = xg - (int64_t)S.t[s].base;
        if (s != 0 && !rv2_a2t(S.t[s], x, &x, &flip)) return false;
        if (tt != 0) { int64_t y; if (!rv2_t2a(S.t[tt], x, &y, &flip)) return false; yg = (int64_t)S.t[tt].base + y; } else yg = x; }
    const uint32_t nd = g.dir ^ (uint32_t)flip; ng->dir = nd;
    ng->c = nd ? (uint64_t)(yg - 1 + (int64_t)o) : (uint64_t)(yg + 1 - (int64_t)o);
#ifdef RV2_FAULT_LIFT
    ng->c += 1;                                                     // negative control: the lifted diagonal one base off
#endif
    return true;
}

// ---------------------------------------------------------------- block decode (the r3_decode_stream of v1 + target ids + LIFT)
// The space is reached through `Acc`: either `ByteAcc` below (S.G, one byte per base, the original path) or
// `AccSpace` of rv2_acc.h (the anchors two bits per base). Only the literal context reads the space here; everything
// else in this function is coordinates.
struct ByteAcc {
    const uint8_t* g = nullptr;
    inline uint8_t at(uint64_t x) const { return g[x]; }
};
RR_HD static inline int r3_refbase_byteacc(const ByteAcc& A, uint64_t ref_n, R3Diag g, uint32_t o) {
    const int64_t q = g.dir ? (int64_t)g.c - (int64_t)o : (int64_t)g.c + (int64_t)o;
    if (q < 0 || (uint64_t)q >= ref_n) return 4;
    const uint8_t b = g.dir ? rr_comp(A.at((uint64_t)q)) : A.at((uint64_t)q);
    const int k = r3_b2(b);
    return k < 4 ? k : 4;
}
static inline int rv2_refbase(const ByteAcc& A, uint64_t ref_n, R3Diag g, uint32_t o) { return r3_refbase_byteacc(A, ref_n, g, o); }

template <class Sink, class Acc>
static inline int r4_decode_stream_acc(const uint8_t* src, uint32_t n, const R4Tab* T, const Space& S, const Acc& A, uint32_t blen, R3Diag start, Sink& sk) {
    const uint64_t ref_n = S.Gn; const int NT = (int)S.t.size();
    R3Dec d; r3_dinit(&d, src, n); if (d.bad) return -1;
    R3Diag cache[4]; int nc = 1; cache[0] = start; int prevk = -1; uint32_t o = 0, k = 0, li = 0;
    const uint8_t B[6] = {'A', 'C', 'G', 'T', 'N', 0};
    while (o < blen) {
        const uint64_t ll = r4_dval(&d, T, R4_LL + r4_kclass(prevk));
        if (d.bad || ll > blen - o) return -1;
        if (ll) {
            uint32_t p1 = 16, p2 = 16;
            for (uint32_t j = 0; j < ll; j++) {
                int c;
                if (ll <= 8) c = R4_LITR + rv2_refbase(A, ref_n, cache[0], o + j) * 2 + (j == 0 ? 1 : 0);
                else c = R4_LIT2 + (p1 == 16 ? 16 : (int)(p2 == 16 ? p1 : p2 * 4 + p1) % 16);
                const uint32_t s = r4_dsym(&d, T, c);
                if (!sk.lit(o + j, li + j, s < 5 ? B[s] : (uint8_t)r3_draw(&d, 8))) return -1;
                p2 = p1; p1 = s < 4 ? s : 0;
            }
            if (!sk.op(0, li, o, (uint32_t)ll)) return -1;
            k++; li += (uint32_t)ll; o += (uint32_t)ll;
        }
        if (o == blen) break;
        const int kind = (int)r4_dsym(&d, T, R4_KIND + r3_llc(ll) * 4 + r4_kclass(prevk));
        uint64_t p = 0; uint32_t dir = 0; int ki = 0; int64_t dd = 0; uint64_t dist = 0; int tt = 0;
        if (kind == K_CONT || kind == K_DELTA) { if (kind == K_DELTA) { const uint64_t z = r4_dval(&d, T, R4_DELTA + (ll == 0 ? 0 : 1)); dd = (int64_t)(z >> 1) ^ -(int64_t)(z & 1); } }
        else if (kind == K_REP) { ki = 1 + (int)r4_dsym(&d, T, R4_REPK); const uint64_t z = r4_dval(&d, T, R4_REPD); dd = (int64_t)(z >> 1) ^ -(int64_t)(z & 1); if (ki >= nc) return -1; }
        else if (kind == K_ABS) { dir = r4_dsym(&d, T, R4_DIR); tt = (int)r4_dsym(&d, T, R4_ABST); p = r3_draw(&d, 32); if (tt >= NT) return -1; p += S.t[tt].base;
#ifdef RV2_FAULT_ABS
            if (tt > 0) p += 1;                                     // negative control: ABS into an anchor one base off
#endif
#ifdef RV2_FAULT_REF
            if (tt == 0) p += 1;                                    // negative control: ABS into the reference one base off
#endif
        }
        else if (kind == K_SELF) { dist = r4_dval(&d, T, R4_SELF);
#ifdef RV2_FAULT_SELF
            dist += 1;                                              // negative control: a self copy resolves one base earlier
#endif
        }
        else if (kind == K_FLIP) { ki = (int)r4_dsym(&d, T, R4_FLIPK); const uint64_t z = r4_dval(&d, T, R4_FLIPD); dd = (int64_t)(z >> 1) ^ -(int64_t)(z & 1); if (ki >= nc) return -1; }
        else if (kind == K_LIFT) { tt = (int)r4_dsym(&d, T, R4_LIFTT); const uint64_t z = r4_dval(&d, T, R4_LIFTD + (ll == 0 ? 0 : 1)); dd = (int64_t)(z >> 1) ^ -(int64_t)(z & 1); if (tt >= NT) return -1; }
        else return -1;
        const uint64_t L = r4_dval(&d, T, R4_LEN + kind) + RR_MINL;
        if (d.bad || L > blen - o) return -1;
        uint32_t okind; uint64_t osrc;
        if (kind == K_SELF) { if (dist == 0 || dist > o) return -1; okind = 2; osrc = o - dist; }
        else {
            if (kind == K_FLIP) { const R3Diag g = cache[ki]; dir = g.dir ^ 1; const int64_t e = r3_locus(g, o) + dd; if (e < 0) return -1; p = (uint64_t)e; }
            else if (kind == K_LIFT) { R3Diag ng; if (!rv2_lift(S, cache[0], o, tt, &ng)) return -1; dir = ng.dir; const int64_t e = r3_exp(ng, o, L) + dd; if (e < 0) return -1; p = (uint64_t)e; }
            else if (kind != K_ABS) { const R3Diag g = cache[kind == K_REP ? ki : 0]; dir = g.dir; const int64_t e = r3_exp(g, o, L) + dd; if (e < 0) return -1; p = (uint64_t)e; }
            if (p > ref_n || L > ref_n - p || !rv2_inside(S, p, L)) return -1;
            okind = dir ? 3 : 1; osrc = p; r3_push(cache, &nc, r3_diag(dir, p, o, L));
        }
        if (!sk.op(okind, osrc, o, (uint32_t)L)) return -1;
        k++; o += (uint32_t)L; prevk = kind;
    }
    if (d.x != R3_L || d.p != d.e) return -1;
    return (int)k;
}
template <class Sink>
static inline int r4_decode_stream(const uint8_t* src, uint32_t n, const R4Tab* T, const Space& S, uint32_t blen, R3Diag start, Sink& sk) {
    ByteAcc A; A.g = S.G;
    return r4_decode_stream_acc(src, n, T, S, A, blen, start, sk);
}
static inline int r4_decode_block(const uint8_t* src, uint32_t n, const R4Tab* T, const Space& S, uint32_t blen, R3Diag start, RrOp* ops, uint32_t maxops, uint8_t* lit, uint32_t litcap) {
    R3ArraySink sk; sk.ops = ops; sk.maxops = maxops; sk.k = 0; sk.lbuf = lit; sk.litcap = litcap;
    return r4_decode_stream(src, n, T, S, blen, start, sk);
}
template <class Acc>
static inline int r4_decode_block_acc(const uint8_t* src, uint32_t n, const R4Tab* T, const Space& S, const Acc& A, uint32_t blen, R3Diag start, RrOp* ops, uint32_t maxops, uint8_t* lit, uint32_t litcap) {
    R3ArraySink sk; sk.ops = ops; sk.maxops = maxops; sk.k = 0; sk.lbuf = lit; sk.litcap = litcap;
    return r4_decode_stream_acc(src, n, T, S, A, blen, start, sk);
}
