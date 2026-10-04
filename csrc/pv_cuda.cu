// pv_cuda.cu - panvram on the GPU: the refrel3 queue kernel over a resident cohort, then the output transforms.
// Queue kernel ported from aceapex research/refrel/refrel3_gpu.cu @ 5b6d5ce (r3q_kernel): one warp per refrel3 block;
// lane 0 decodes the block's rANS stream (r3_decode_stream of refrel3.h, the code the CPU path runs), literals straight
// into the block buffer in shared memory and copies into a ring of QN ops; lanes 1..31 take the ops in order as they come
// (reference / reverse complement / self); then the warp writes the window's part. Changes: the block size Q is the
// cohort's (runtime, shared memory Q + ring), many assemblies resident at once (pooled arrays indexed per assembly), a
// window is (assembly, start); every wait is bounded (status 4 instead of a hang).
// Transforms after it: lower-case runs (as in the FASTA), reverse complement (per window), tokens A/a 0 C/c 1 G/g 2 T/t 3 else 4.
#include <torch/extension.h>
#include <c10/cuda/CUDAStream.h>
#include <c10/cuda/CUDAException.h>
#include <c10/cuda/CUDAGuard.h>
#include "refrel3.h"

namespace {
struct QOp { uint32_t src; uint16_t dst, len; uint32_t kind; };
#define QN 64
#define QSPIN (1u << 26)
struct QSink { uint8_t* blk; QOp* q; volatile int* head; volatile int* tail; int* status;
    __device__ bool lit(uint32_t o, uint32_t, uint8_t b) { blk[o] = b; return true; }
    __device__ bool op(uint32_t kind, uint64_t src, uint32_t dst, uint32_t len) {
        if (kind == 0) return true;                                          // literals are already in the block buffer
        const int h = *head; uint32_t spin = 0;
        while (h - *tail >= QN) { if (++spin > QSPIN) { atomicOr(status, 4); return false; } }
        QOp& e = q[h % QN]; e.src = (uint32_t)src; e.dst = (uint16_t)dst; e.len = (uint16_t)len; e.kind = kind;
        __threadfence_block(); *head = h + 1; return true; } };
struct Cohort { const uint8_t* P; const int64_t* pay_base; const uint32_t* off; const int64_t* st; const R3Tab* T; const int64_t* blk_base; const int64_t* nbases;
                const uint8_t* ref; uint64_t ref_n; uint32_t Q; int64_t nasm; };

// status bits: 1 decode error, 4 bounded wait expired, 8 window outside the assembly
__global__ void __launch_bounds__(32, 16) pv_queue_kernel(Cohort D, const int64_t* __restrict__ wasm, const int64_t* __restrict__ wstart, uint32_t W, uint32_t bpw,
                                                          uint8_t* __restrict__ out, int* status) {
    extern __shared__ __align__(16) uint8_t sm[];
    const uint32_t Q = D.Q; uint8_t* blk = sm; QOp* q = (QOp*)(sm + Q);
    __shared__ int s_head, s_tail, s_done, s_nk;
    volatile int* head = &s_head; volatile int* tail = &s_tail; volatile int* done = &s_done;
    const uint32_t lane = threadIdx.x;
    const uint64_t w = blockIdx.x / bpw; const uint32_t j = blockIdx.x % bpw;
    const int64_t a = wasm[w]; if (a < 0 || a >= D.nasm) { if (lane == 0 && j == 0) atomicOr(status, 8); return; }
    const uint64_t s = (uint64_t)wstart[w], n_bases = (uint64_t)D.nbases[a];
    if (wstart[w] < 0 || s + W > n_bases) { if (lane == 0 && j == 0) atomicOr(status, 8); return; }
    const uint64_t b0 = s / Q, b1 = (s + W - 1) / Q, b = b0 + j; if (b > b1) return;
    const uint64_t bst = b * Q; const uint32_t blen = (uint32_t)min((uint64_t)Q, n_bases - bst);
    const uint64_t bb = (uint64_t)D.blk_base[a] + b; const uint8_t* src = D.P + D.pay_base[a] + D.off[bb]; const uint32_t slen = D.off[bb + 1] - D.off[bb];
    R3Diag g0; g0.c = (uint64_t)(D.st[bb] >> 1); g0.dir = (uint32_t)(D.st[bb] & 1);
    if (lane == 0) { *head = 0; *tail = 0; *done = 0; }
    __syncwarp();
    if (lane == 0) {
        QSink sk; sk.blk = blk; sk.q = q; sk.head = head; sk.tail = tail; sk.status = status;
        s_nk = r3_decode_stream(src, slen, D.T + a, D.ref, D.ref_n, blen, g0, sk);
        __threadfence_block(); *done = s_nk < 0 ? 2 : 1;
    } else {
        const unsigned m = 0xFFFFFFFEu; int next = 0; uint32_t spin = 0;
        for (;;) {
            int h = 0, d = 0; if (lane == 1) { d = *done; __threadfence_block(); h = *head; }   // done first: then head is final
            h = __shfl_sync(m, h, 1); d = __shfl_sync(m, d, 1);
            if (next < h) {
                const QOp e = q[next % QN];
                if (e.kind == 1) for (uint32_t i = lane - 1; i < e.len; i += 31) blk[e.dst + i] = D.ref[e.src + i];
                else if (e.kind == 3) for (uint32_t i = lane - 1; i < e.len; i += 31) blk[e.dst + i] = rr_comp(D.ref[e.src + e.len - 1 - i]);
                else { const uint32_t dist = (uint32_t)e.dst - e.src;
                    for (uint32_t c = 0; c < e.len;) { const uint32_t step = min(dist, (uint32_t)e.len - c);
                        for (uint32_t i = lane - 1; i < step; i += 31) blk[e.dst + c + i] = blk[e.src + c + i];
                        __syncwarp(m); c += step; } }
                __syncwarp(m); next++; if (lane == 1) { __threadfence_block(); *tail = next; } spin = 0;
            } else if (d) { if (next >= h) break; }
            else if (++spin > QSPIN) { if (lane == 1) atomicOr(status, 4); break; }
        }
    }
    __syncwarp();
    if (s_nk < 0) { if (lane == 0) atomicOr(status, 1); return; }
    const uint64_t x0 = max(s, bst), x1 = min(s + W, bst + blen);
    for (uint64_t x = x0 + lane; x < x1; x += 32) out[w * W + (x - s)] = blk[x - bst];
}

#define CH 8192u                                                             // positions per CUDA block in the transforms
// lower-case runs: grid (chunks, windows)
__global__ void pv_case_kernel(const int64_t* __restrict__ wasm, const int64_t* __restrict__ wstart, uint32_t W, const int64_t* __restrict__ low_s, const int64_t* __restrict__ low_l,
                               const int64_t* __restrict__ low_base, const int64_t* __restrict__ low_cnt, uint8_t* __restrict__ out) {
    const uint64_t w = blockIdx.y; const int64_t a = wasm[w]; const int64_t cnt = low_cnt[a]; if (!cnt) return;
    const int64_t s = wstart[w], lo = s + (int64_t)blockIdx.x * CH, hi = min(s + (int64_t)W, lo + (int64_t)CH); if (lo >= hi) return;
    const int64_t* S = low_s + low_base[a]; const int64_t* L = low_l + low_base[a];
    __shared__ int64_t first;
    if (threadIdx.x == 0) { int64_t l = 0, h = cnt; while (l < h) { const int64_t m = (l + h) / 2; if (S[m] + L[m] <= lo) l = m + 1; else h = m; } first = l; }
    __syncthreads();
    uint8_t* o = out + w * W;
    for (int64_t r = first; r < cnt && S[r] < hi; r++) { const int64_t x0 = max(lo, S[r]), x1 = min(hi, S[r] + L[r]);
        for (int64_t x = x0 + threadIdx.x; x < x1; x += blockDim.x) o[x - s] |= 0x20; }
}
__device__ __forceinline__ uint8_t comp_case(uint8_t c) { switch (c) { case 'A': return 'T'; case 'C': return 'G'; case 'G': return 'C'; case 'T': return 'A';
    case 'a': return 't'; case 'c': return 'g'; case 'g': return 'c'; case 't': return 'a'; default: return c; } }
__device__ __forceinline__ uint8_t token(uint8_t c) { switch (c | 0x20) { case 'a': return 0; case 'c': return 1; case 'g': return 2; case 't': return 3; default: return 4; } }
// reverse complement (rc[w]) and / or tokens over pairs (i, W - 1 - i): grid (chunks of the first half, windows)
__global__ void pv_rc_tok_kernel(uint32_t W, const uint8_t* __restrict__ rc, bool tokens, uint8_t* __restrict__ out) {
    const uint64_t w = blockIdx.y; const bool r = rc && rc[w]; if (!r && !tokens) return;
    uint8_t* o = out + w * W; const uint32_t half = (W + 1) / 2, i0 = blockIdx.x * CH, i1 = min(half, i0 + CH);
    for (uint32_t i = i0 + threadIdx.x; i < i1; i += blockDim.x) { const uint32_t j = W - 1 - i; uint8_t x = o[i], y = o[j];
        if (r) { const uint8_t t = comp_case(x); x = comp_case(y); y = t; }
        if (tokens) { x = token(x); y = token(y); }
        o[i] = x; if (j != i) o[j] = y; }
}
}  // namespace

std::vector<torch::Tensor> pv_windows_cuda(torch::Tensor P, torch::Tensor pay_base, torch::Tensor off, torch::Tensor st, torch::Tensor tab, torch::Tensor blk_base,
                                           torch::Tensor nbases, torch::Tensor low_s, torch::Tensor low_l, torch::Tensor low_base, torch::Tensor low_cnt, torch::Tensor ref,
                                           int64_t Q, torch::Tensor wasm, torch::Tensor wstart, int64_t W, c10::optional<torch::Tensor> rc, bool tokens, bool apply_case) {
    TORCH_CHECK(P.is_cuda(), "pooled cohort not on a CUDA device");
    TORCH_CHECK(Q == 1024 || Q == 2048 || Q == 4096 || Q == 16384, "block size");
    TORCH_CHECK(W >= 0 && W < (1ll << 31), "W");
    const c10::cuda::CUDAGuard guard(P.device());
    auto dev = P.device();
    wasm = wasm.to(dev, torch::kInt64).contiguous(); wstart = wstart.to(dev, torch::kInt64).contiguous();
    const int64_t n = wasm.numel(); TORCH_CHECK(wstart.numel() == n, "assembly and start counts differ");
    torch::Tensor out = torch::empty({n, W}, torch::TensorOptions().dtype(torch::kUInt8).device(dev));
    torch::Tensor status = torch::zeros({1}, torch::TensorOptions().dtype(torch::kInt32).device(dev));
    if (n == 0 || W == 0) return {out, status};
    torch::Tensor r; if (rc && rc->defined()) { r = rc->to(dev, torch::kUInt8).contiguous(); TORCH_CHECK(r.numel() == n, "reverse_complement mask size"); }
    auto stream = c10::cuda::getCurrentCUDAStream();
    Cohort D; D.P = P.data_ptr<uint8_t>(); D.pay_base = pay_base.data_ptr<int64_t>(); D.off = (const uint32_t*)off.data_ptr<int32_t>(); D.st = st.data_ptr<int64_t>();
    D.T = (const R3Tab*)tab.data_ptr<uint8_t>(); D.blk_base = blk_base.data_ptr<int64_t>(); D.nbases = nbases.data_ptr<int64_t>();
    D.ref = ref.data_ptr<uint8_t>(); D.ref_n = (uint64_t)ref.numel(); D.Q = (uint32_t)Q; D.nasm = nbases.numel();
    const uint64_t bpw = ((uint64_t)W + Q - 1) / Q + 1; TORCH_CHECK((uint64_t)n * bpw < (1ull << 31), "too many blocks in one call (windows x blocks per window >= 2^31)");
    const size_t smem = (size_t)Q + QN * sizeof(QOp);
    pv_queue_kernel<<<(unsigned)(n * bpw), 32, smem, stream>>>(D, wasm.data_ptr<int64_t>(), wstart.data_ptr<int64_t>(), (uint32_t)W, (uint32_t)bpw, out.data_ptr<uint8_t>(), status.data_ptr<int>());
    C10_CUDA_KERNEL_LAUNCH_CHECK();
    // grid.y is limited to 65535: transforms in slices of windows
    for (int64_t w0 = 0; w0 < n; w0 += 65535) { const int64_t k = std::min<int64_t>(65535, n - w0);
        if (apply_case && !tokens && low_s.numel()) {
            dim3 g((unsigned)((W + CH - 1) / CH), (unsigned)k);
            pv_case_kernel<<<g, 256, 0, stream>>>(wasm.data_ptr<int64_t>() + w0, wstart.data_ptr<int64_t>() + w0, (uint32_t)W, low_s.data_ptr<int64_t>(), low_l.data_ptr<int64_t>(),
                                                  low_base.data_ptr<int64_t>(), low_cnt.data_ptr<int64_t>(), out.data_ptr<uint8_t>() + w0 * W);
            C10_CUDA_KERNEL_LAUNCH_CHECK(); }
        if (r.defined() || tokens) {
            dim3 g((unsigned)(((W + 1) / 2 + CH - 1) / CH), (unsigned)k);
            pv_rc_tok_kernel<<<g, 256, 0, stream>>>((uint32_t)W, r.defined() ? r.data_ptr<uint8_t>() + w0 : nullptr, tokens, out.data_ptr<uint8_t>() + w0 * W);
            C10_CUDA_KERNEL_LAUNCH_CHECK(); }
    }
    return {out, status};
}
