import re
from pathlib import Path
cu = Path('/content/panvram/csrc/pv_cuda.cu'); core = Path('/content/panvram/csrc/pv_core.cpp')
s, c = cu.read_text(), core.read_text()

A1 = 'const uint32_t lane = threadIdx.x;'
A2 = '        s_nk = r3_decode_stream(src, slen, D.T + a, ref, D.ref_n, blen, g0, sk);'
A3 = '        const unsigned m = 0xFFFFFFFEu; int next = 0; uint32_t spin = 0;'
A4 = '''                __syncwarp(m); next++; if (lane == 1) { __threadfence_block(); *tail = next; } spin = 0;'''
A5 = '            if (next < h) {'
A6 = '            } else if (d) { if (next >= h) break; }'
for i, a in enumerate((A1, A2, A3, A4, A5, A6), 1):
    assert s.count(a) == 1, f'anchor A{i} found {s.count(a)} times - upstream changed, stopping'

A0 = 'template <bool PACKED, bool COMPACT>'
assert s.count(A0) == 1, f'anchor A0 found {s.count(A0)} times - upstream changed, stopping'
s = s.replace(A0, '__device__ unsigned long long pv_cyc[4];\n\n' + A0, 1)
assert 'pv_cyc[4]' in s
s = s.replace(A2, '        const unsigned long long _t0 = clock64();\n' + A2 +
              '\n        atomicAdd(&pv_cyc[0], clock64() - _t0); atomicAdd(&pv_cyc[3], 1ull);', 1)
s = s.replace(A3, A3 + '\n        unsigned long long _cwork = 0; const unsigned long long _cl0 = clock64();', 1)
s = s.replace(A5, A5 + '\n                const unsigned long long _w0 = clock64();', 1)
s = s.replace(A4, '                if (lane == 1) _cwork += clock64() - _w0;\n' + A4, 1)
s = s.replace(A6, A6.replace('break; }', 'break; }') +
              '\n            // accounting below on leaving the loop', 1)
s = s.replace('''        }
    }
    __syncwarp();
    if (s_nk < 0)''', '''        }
        if (lane == 1) { atomicAdd(&pv_cyc[1], _cwork); atomicAdd(&pv_cyc[2], clock64() - _cl0); }
    }
    __syncwarp();
    if (s_nk < 0)''', 1)

READER = '''
#include <vector>
std::vector<unsigned long long> pv_read_cycles(bool reset) {
    unsigned long long h[4] = {0, 0, 0, 0};
    cudaMemcpyFromSymbol(h, pv_cyc, sizeof(h));
    if (reset) { unsigned long long z[4] = {0, 0, 0, 0}; cudaMemcpyToSymbol(pv_cyc, z, sizeof(z)); }
    return std::vector<unsigned long long>(h, h + 4);
}
'''
s = s + READER
assert 'PYBIND11_MODULE' in c or 'm.def(' in c, 'pv_core.cpp does not look like the pybind module'
hook = re.search(r'(PYBIND11_MODULE\s*\([^)]*\)\s*\{)', c)
assert hook, 'could not find the pybind module block'
c = (c[:hook.end()] +
     '\n#ifdef PV_WITH_CUDA\n    { extern std::vector<unsigned long long> pv_read_cycles(bool);'
     '\n      m.def("read_cycles", &pv_read_cycles); }\n#endif\n' +
     c[hook.end():])
cu.write_text(s); core.write_text(c)
print('patched; diff summary:')
!cd /content/panvram && git diff --stat -- csrc/
!cd /content/panvram && git diff -- csrc/pv_cuda.cu | head -60