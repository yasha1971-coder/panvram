// rv2_acc.h - a copy of refrel3.h's reference access on an accessor, so the eight anchor assemblies of the refrel v2
// prototype can be held two bits per base instead of one byte.
//
// What is copied, and from where, unchanged except that the reference is reached through `A` instead of a byte
// pointer: `r3_refbase` of tools/rr3bench/refrel3.h and `exec_fast` of tools/rr3bench/refrel.cpp. The shared v1 code
// itself is NOT modified: tools/rr3bench and tools/window_cost keep the binaries their logs record.
//
// Layout of the target space. Target 0 is the reference (T2T), still one byte per base: it is the source of every v1
// archive and of the anchors' own copies, and it carries no non-ACGT base at all. Targets 1..K are the anchors,
// packed: 16 bases per uint32 word, A 0 C 1 G 2 T 3, plus a sorted list of runs of any other byte (start, length,
// byte) - 17 bytes per run, 95 runs over the eight anchors of the N = 558 build.
//
// Reads are by range, not by base: `copy` and `copy_rc` unpack whole words and then patch the few exception runs that
// overlap the range, one lower_bound per range. A per-base binary search would cost more than the packing saves.
#pragma once
#include <algorithm>
#include <cstring>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

#include "refrel3.h"

// The global coordinate system is the one the encoder used: target 0 at [0, rn), then the anchors back to back with
// no gap. Packed storage is separate and 16-aligned per anchor, so a block never shares a word with another block:
// `ps[k]` is where anchor k starts inside the word array. Mixing the two was a real bug, caught by the judge.
struct AccSpace {
    const uint8_t* ref = nullptr;                 // target 0, bytes
    uint64_t rn = 0;
    const uint32_t* w = nullptr;                  // packed anchors
    uint64_t wn = 0;
    int nt = 0;                                   // anchors
    uint64_t gb[33] = {0};                        // anchor k: global base, length, packed start
    uint64_t gl[33] = {0};
    uint64_t ps[33] = {0};
    const uint64_t* es = nullptr;                 // exception runs in GLOBAL coordinates, sorted, disjoint
    const uint64_t* el = nullptr;
    const uint8_t* eb = nullptr;
    int64_t ne = 0;

    static inline uint8_t unpack(uint32_t word, uint32_t slot) {
        static const char B[4] = {'A', 'C', 'G', 'T'};
        return (uint8_t)B[(word >> (slot * 2)) & 3u];
    }

    inline int target_of(uint64_t x) const {      // which anchor holds global x; -1 if it is in target 0
        for (int k = nt - 1; k >= 0; k--) if (x >= gb[k]) return x < gb[k] + gl[k] ? k : -2;
        return -1;
    }

    inline int64_t exc_lo(uint64_t x) const {     // first exception run that can touch x
        int64_t lo = 0, hi = ne;
        while (lo < hi) { const int64_t m = (lo + hi) / 2; if (es[m] + el[m] <= x) lo = m + 1; else hi = m; }
        return lo;
    }

    inline uint8_t at(uint64_t x) const {
        if (x < rn) return ref[x];
        if (ne) { const int64_t i = exc_lo(x); if (i < ne && es[i] <= x) return eb[i]; }
        const int k = target_of(x);
        if (k < 0) return 'N';
        const uint64_t p = ps[k] + (x - gb[k]);
        return unpack(w[p >> 4], (uint32_t)(p & 15));
    }

    // out[0, len) = the bases at [x, x + len); the range never crosses a target (rv2_inside guarantees it)
    inline void copy(uint64_t x, uint32_t len, uint8_t* out) const {
        if (x + len <= rn) { memcpy(out, ref + x, len); return; }
        const int k = target_of(x);
        if (k < 0) { memset(out, 'N', len); return; }
        uint64_t p = ps[k] + (x - gb[k]);
        uint32_t i = 0;
        while (i < len) {
            const uint32_t slot = (uint32_t)(p & 15);
            const uint32_t take = std::min<uint32_t>(16u - slot, len - i);
            const uint32_t word = w[p >> 4];
            for (uint32_t j = 0; j < take; j++) out[i + j] = unpack(word, slot + j);
            i += take; p += take;
        }
        if (!ne) return;
        for (int64_t e = exc_lo(x); e < ne && es[e] < x + len; e++) {
            const uint64_t a = std::max(es[e], x), b = std::min(es[e] + el[e], x + len);
            if (a < b) memset(out + (a - x), eb[e], (size_t)(b - a));
        }
    }

    // out[j] = comp(base at x + len - 1 - j)
    inline void copy_rc(uint64_t x, uint32_t len, uint8_t* out) const {
        if (x + len <= rn) {
            const uint8_t* r = ref + x + len - 1;
            for (uint32_t j = 0; j < len; j++) out[j] = COMP[*(r - j)];
            return;
        }
        copy(x, len, out);
        for (uint32_t j = 0; j < len / 2; j++) {
            const uint8_t a = COMP[out[j]], b = COMP[out[len - 1 - j]];
            out[j] = b; out[len - 1 - j] = a;
        }
        if (len & 1) out[len / 2] = COMP[out[len / 2]];
    }
};

// = r3_refbase of tools/rr3bench/refrel3.h, with the reference reached through A
RR_HD static inline int r3_refbase_acc(const AccSpace& A, uint64_t ref_n, R3Diag g, uint32_t o) {
    const int64_t q = g.dir ? (int64_t)g.c - (int64_t)o : (int64_t)g.c + (int64_t)o;
    if (q < 0 || (uint64_t)q >= ref_n) return 4;
    const uint8_t b = g.dir ? rr_comp(A.at((uint64_t)q)) : A.at((uint64_t)q);
    const int k = r3_b2(b);
    return k < 4 ? k : 4;
}

static inline int rv2_refbase(const AccSpace& A, uint64_t ref_n, R3Diag g, uint32_t o) { return r3_refbase_acc(A, ref_n, g, o); }

// = exec_fast of tools/rr3bench/refrel.cpp, with the reference reached through A
static inline void exec_fast_acc(const RrOp* ops, int n, const uint8_t* lit, const AccSpace& A, uint8_t* out) {
    for (int k = 0; k < n; k++) {
        const RrOp& q = ops[k];
        if (q.kind == 0) memcpy(out + q.dst, lit + q.src, q.len);
        else if (q.kind == 1) A.copy(q.src, q.len, out + q.dst);
        else if (q.kind == 2) {
            const uint32_t dist = q.dst - (uint32_t)q.src;
            if (dist >= q.len) memcpy(out + q.dst, out + q.src, q.len);
            else for (uint32_t j = 0; j < q.len; j++) out[q.dst + j] = out[q.src + j];
        } else A.copy_rc(q.src, q.len, out + q.dst);
    }
}

// Builder. Q is a multiple of 16, so two blocks never share a word and threads can pack in parallel without a lock.
// `acc_pack_block` writes the 2-bit codes and reports the runs of any other byte, in packed-region coordinates; the
// caller turns them into global coordinates and merges them into one sorted list.
static inline void acc_pack_block(uint32_t* w, uint64_t p0, const uint8_t* bases, uint32_t n,
                                  std::vector<uint64_t>& run_start, std::vector<uint64_t>& run_len,
                                  std::vector<uint8_t>& run_byte) {
    // Every block must start on a word boundary, or two threads would read-modify-write the same word: the caller
    // pads each anchor's start to a multiple of 16 and Q is a multiple of 16, so only an anchor's last block ends
    // mid-word and no other thread touches that word.
    if (p0 % 16) { fprintf(stderr, "acc_pack_block: block start %llu is not a multiple of 16\n", (unsigned long long)p0); abort(); }
    uint32_t j = 0;
    while (j < n) {
        const uint64_t p = p0 + j;
        const uint32_t slot = (uint32_t)(p & 15);
        const uint32_t take = std::min<uint32_t>(16u - slot, n - j);
        uint32_t word = 0;
        for (uint32_t k = 0; k < take; k++) {
            const uint8_t b = bases[j + k];
            const uint32_t code = b == 'A' ? 0u : b == 'C' ? 1u : b == 'G' ? 2u : b == 'T' ? 3u : 0u;
            word |= code << ((slot + k) * 2);
        }
        if (slot == 0 && take == 16) w[p >> 4] = word;             // whole word: plain store
        else w[p >> 4] |= word;                                     // partial word at a block edge: Q % 16 == 0 keeps this race-free
        j += take;
    }
    for (uint32_t i = 0; i < n;) {                                  // runs of anything that is not A/C/G/T
        const uint8_t b = bases[i];
        if (b == 'A' || b == 'C' || b == 'G' || b == 'T') { i++; continue; }
        uint32_t e = i + 1;
        while (e < n && bases[e] == b) e++;
        run_start.push_back(p0 + i); run_len.push_back(e - i); run_byte.push_back(b);
        i = e;
    }
}

// Merge per-thread run lists into one sorted, disjoint list (sorted by start; the lists never overlap).
static inline void acc_merge_runs(std::vector<uint64_t>& es, std::vector<uint64_t>& el, std::vector<uint8_t>& eb) {
    const size_t n = es.size();
    std::vector<size_t> idx(n);
    for (size_t i = 0; i < n; i++) idx[i] = i;
    std::sort(idx.begin(), idx.end(), [&](size_t a, size_t b) { return es[a] < es[b]; });
    std::vector<uint64_t> s2(n), l2(n); std::vector<uint8_t> b2(n);
    for (size_t i = 0; i < n; i++) { s2[i] = es[idx[i]]; l2[i] = el[idx[i]]; b2[i] = eb[idx[i]]; }
    es.swap(s2); el.swap(l2); eb.swap(b2);
}
