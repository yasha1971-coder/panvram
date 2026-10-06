// pv_compact.h - compact resident forms (opt-in; the archives and their format are unchanged), shared by the CPU path
// and the CUDA kernel:
//   PackedRef   the decoded reference in 2 bits per base (16 bases per uint32 word, A 0 C 1 G 2 T 3) + a sorted table
//               of runs of any other byte (N, IUPAC, ...): (start, length, byte); an index returns the same byte the
//               plain byte array would
//   CompactTab  the block table (payload offset + start state per block) in two levels: per group of PV_G blocks the
//               offset and packed state of its first block; per block a 16-bit payload length and a 16-bit code
//               (state delta from the chain prediction c +- Q, x 2 + strand; PV_EXC = look the state up in a sorted
//               exception table). Offset and state of block b: one pass over at most PV_G - 1 earlier blocks of its group.
#pragma once
#include <stdint.h>
#include "refrel3.h"

namespace pv {

#define PV_G 32
#define PV_EXC ((int16_t)-32768)

struct PackedRef {
    const uint32_t* w; uint64_t n;
    const int64_t* es; const int64_t* el; const uint8_t* eb; int64_t ne;     // exception runs, sorted, disjoint
    RR_HD uint8_t operator[](uint64_t q) const {
        if (ne) { int64_t lo = 0, hi = ne; while (lo < hi) { const int64_t m = (lo + hi) / 2; if ((uint64_t)(es[m] + el[m]) <= q) lo = m + 1; else hi = m; }
            if (lo < ne && (uint64_t)es[lo] <= q) return eb[lo]; }
        const uint32_t v = (w[q >> 4] >> ((q & 15) * 2)) & 3;
        return (uint8_t)(v == 0 ? 'A' : v == 1 ? 'C' : v == 2 ? 'G' : 'T');
    }
};

struct CompactTab {
    const int64_t* grp_off; const int64_t* grp_st;   // per group (pooled): payload offset of its first block (relative to the assembly), packed state
    const uint16_t* len; const int16_t* code;        // per block (pooled)
    const int64_t* exc_b; const int64_t* exc_st; int64_t nexc;   // pooled block index -> packed state, sorted
    uint32_t Q;
    // block b of the assembly whose first pooled block is gb0 and first group g0
    RR_HD void get(uint64_t gb0, uint64_t g0, uint64_t b, uint32_t* off, uint32_t* blen, R3Diag* st) const {
        const uint64_t g = g0 + b / PV_G, first = b - b % PV_G;
        int64_t o = grp_off[g]; int64_t s = grp_st[g];
        for (uint64_t k = first + 1; k <= b; k++) {
            o += len[gb0 + k - 1];
            const int16_t c = code[gb0 + k];
            if (c == PV_EXC) { int64_t lo = 0, hi = nexc; const int64_t key = (int64_t)(gb0 + k); while (lo < hi) { const int64_t m = (lo + hi) / 2; if (exc_b[m] < key) lo = m + 1; else hi = m; } s = exc_st[lo]; }
            else { const int64_t pc = s >> 1; const int64_t pd = s & 1; const int64_t pred = pd ? pc - (int64_t)Q : pc + (int64_t)Q;
                const int64_t nc = pred + (c >> 1); s = (int64_t)((uint64_t)nc << 1) | (c & 1); }
        }
        *off = (uint32_t)o; *blen = len[gb0 + b]; st->c = (uint64_t)(s >> 1); st->dir = (uint32_t)(s & 1);
    }
};

}  // namespace pv
