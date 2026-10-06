// pv_core.cpp - panvram host side and the Python binding (module panvram._C).
// Core: the reference (decoded, upper case, SHA-256) and the cohort's refrel3 v1 archives opened with every check
// (pv_v1.h), their arrays pooled into tensors - payloads, block offsets (32-bit, per assembly), block start states,
// model tables, lower-case runs, contig tables - which Python moves to the device as they are; the CPU decoder of
// windows (the reference path of the CUDA queue kernel, pv_cuda.cu), full decode -> FASTA with all hashes checked.
#include <torch/extension.h>
#include "pv_v1.h"
#include <atomic>
#include <mutex>
#include <thread>

namespace pv {

template <class F> static void parallel(uint64_t n, int T, F&& f) {
    T = (int)std::max<uint64_t>(1, std::min<uint64_t>((uint64_t)std::max(1, T), n)); if (!n) return;
    std::atomic<uint64_t> next{0}; std::vector<std::thread> th;
    for (int t = 0; t < T; t++) th.emplace_back([&] { for (;;) { const uint64_t i = next.fetch_add(1); if (i >= n) break; f(i); } });
    for (auto& x : th) x.join();
}
static int default_threads() { const unsigned h = std::thread::hardware_concurrency(); return h ? (int)h : 4; }

struct Core {
    std::string ref_path; uint8_t ref_sha[32]; uint64_t ref_n = 0; torch::Tensor ref;          // ref: uint8 [ref_n]
    std::vector<std::unique_ptr<Archive>> A; uint32_t Q = 0; int threads = 0;
    torch::Tensor P, pay_base, off, st, tab, blk_base, nbases, low_s, low_l, low_base, low_cnt;  // pooled (CPU)
    torch::Tensor ctg_asm, ctg_local, ctg_boff, ctg_len;                                         // all contigs, cohort order
    // opt-in compact forms (pv_compact.h); the byte reference / full block table are then not kept
    bool packed = false, compact = false; torch::Tensor pw, pes, pel, peb; PackedRef pr{};
    torch::Tensor grp_off, grp_st, blen, bcode, exc_b, exc_st, grp_base, blk0; CompactTab ct{};

    Core(const std::string& reference, const std::vector<std::string>& paths, int T, bool packed_ref, bool compact_blocks)
        : ref_path(reference), threads(T > 0 ? T : default_threads()), packed(packed_ref), compact(compact_blocks) {
        if (paths.empty()) throw Err("no archives");
        { FILE* f = fopen(reference.c_str(), "rb"); if (!f) throw Err(reference + ": cannot open");
          fseek(f, 0, SEEK_END); const long s = ftell(f); fseek(f, 0, SEEK_SET);
          torch::Tensor buf = torch::empty({(int64_t)s}, torch::kUInt8); const size_t r = fread(buf.data_ptr<uint8_t>(), 1, (size_t)s, f); fclose(f);
          if (r != (size_t)s) throw Err(reference + ": short read");
          ref_n = fasta_bases_inplace(buf.data_ptr<uint8_t>(), (uint64_t)s, reference); ref = buf.narrow(0, 0, (int64_t)ref_n);
          rr_sha256(ref.data_ptr<uint8_t>(), ref_n, ref_sha); }
        if (ref_n >= (1ull << 32)) throw Err("reference of 2^32 bases or more: the queue kernel's copy positions are 32-bit");
        if (packed) pack_reference();
        A.resize(paths.size()); std::mutex mu; std::string first;
        parallel(paths.size(), threads, [&](uint64_t i) {
            try { const std::vector<uint8_t> f = slurp(paths[i]); std::unique_ptr<Archive> X(new Archive()); open_v1(f, paths[i], ref_sha, ref_n, *X); A[i] = std::move(X); }
            catch (const std::exception& e) { std::lock_guard<std::mutex> g(mu); if (first.empty()) first = e.what(); } });
        if (!first.empty()) throw Err(first);
        Q = A[0]->Q; for (auto& x : A) if (x->Q != Q) throw Err(x->path + ": block size " + std::to_string(x->Q) + " differs from " + std::to_string(Q) + " (one dataset per cohort)");
        pool();
    }
    // 2 bits per base + runs of every other byte; the byte array is released
    void pack_reference() {
        const uint8_t* R = ref.data_ptr<uint8_t>(); const int64_t nw = (int64_t)((ref_n + 15) / 16);
        pw = torch::zeros({nw}, torch::kInt32); uint32_t* w = (uint32_t*)pw.data_ptr<int32_t>();
        std::vector<int64_t> es, el; std::vector<uint8_t> eb;
        for (uint64_t q = 0; q < ref_n;) { const uint8_t c = R[q]; const int v = c == 'A' ? 0 : c == 'C' ? 1 : c == 'G' ? 2 : c == 'T' ? 3 : -1;
            if (v >= 0) { w[q >> 4] |= (uint32_t)v << ((q & 15) * 2); q++; continue; }
            uint64_t e = q; while (e < ref_n && R[e] == c) e++; es.push_back((int64_t)q); el.push_back((int64_t)(e - q)); eb.push_back(c); q = e; }
        pes = torch::tensor(es.empty() ? std::vector<int64_t>{0} : es, torch::kInt64); pel = torch::tensor(el.empty() ? std::vector<int64_t>{0} : el, torch::kInt64);
        peb = torch::empty({(int64_t)std::max<size_t>(1, eb.size())}, torch::kUInt8); if (!eb.empty()) memcpy(peb.data_ptr<uint8_t>(), eb.data(), eb.size());
        pr.w = w; pr.n = ref_n; pr.es = pes.data_ptr<int64_t>(); pr.el = pel.data_ptr<int64_t>(); pr.eb = peb.data_ptr<uint8_t>(); pr.ne = (int64_t)es.size();
        ref = torch::Tensor();
    }
    template <class F> auto with_ref(F&& f) const { if (packed) return f(pr); const uint8_t* R = ref.data_ptr<uint8_t>(); return f(R); }
    void pool() {
        const int64_t n = (int64_t)A.size(); uint64_t sp = 0, sb = 0, sl = 0, sc = 0;
        std::vector<uint64_t> pb(n), bb(n), lb(n), cb(n);
        for (int64_t i = 0; i < n; i++) { pb[i] = sp; bb[i] = sb; lb[i] = sl; cb[i] = sc; sp += A[i]->payload_bytes; sb += A[i]->nb + 1; sl += A[i]->low_s.size(); sc += A[i]->rec.size(); }
        auto i64 = torch::kInt64;
        P = torch::empty({(int64_t)sp + 64}, torch::kUInt8); P.narrow(0, (int64_t)sp, 64).zero_();
        if (!compact) { off = torch::empty({(int64_t)sb}, torch::kInt32); st = torch::empty({(int64_t)sb}, i64); }
        tab = torch::empty({n * (int64_t)sizeof(R3Tab)}, torch::kUInt8);
        pay_base = torch::empty({n}, i64); blk_base = torch::empty({n}, i64); nbases = torch::empty({n}, i64); low_base = torch::empty({n}, i64); low_cnt = torch::empty({n}, i64);
        low_s = torch::empty({(int64_t)sl}, i64); low_l = torch::empty({(int64_t)sl}, i64);
        ctg_asm = torch::empty({(int64_t)sc}, i64); ctg_local = torch::empty({(int64_t)sc}, i64); ctg_boff = torch::empty({(int64_t)sc}, i64); ctg_len = torch::empty({(int64_t)sc}, i64);
        uint8_t* Pp = P.data_ptr<uint8_t>(); uint32_t* offp = compact ? nullptr : (uint32_t*)off.data_ptr<int32_t>(); int64_t* stp = compact ? nullptr : st.data_ptr<int64_t>(); R3Tab* tp = (R3Tab*)tab.data_ptr<uint8_t>();
        if (compact) build_compact();
        parallel((uint64_t)n, threads, [&](uint64_t i) {
            Archive& X = *A[i];
            memcpy(Pp + pb[i], X.P_own.data(), X.payload_bytes);
            if (!compact) { memcpy(offp + bb[i], X.off_own.data(), (X.nb + 1) * 4); memcpy(stp + bb[i], X.st_own.data(), (X.nb + 1) * 8); }
            memcpy(tp + i, X.T_own.get(), sizeof(R3Tab));
            pay_base.data_ptr<int64_t>()[i] = (int64_t)pb[i]; blk_base.data_ptr<int64_t>()[i] = (int64_t)bb[i]; nbases.data_ptr<int64_t>()[i] = (int64_t)X.nbases;
            low_base.data_ptr<int64_t>()[i] = (int64_t)lb[i]; low_cnt.data_ptr<int64_t>()[i] = (int64_t)X.low_s.size();
            for (size_t r = 0; r < X.low_s.size(); r++) { low_s.data_ptr<int64_t>()[lb[i] + r] = (int64_t)X.low_s[r]; low_l.data_ptr<int64_t>()[lb[i] + r] = (int64_t)X.low_l[r]; }
            for (size_t r = 0; r < X.rec.size(); r++) { const uint64_t k = cb[i] + r; ctg_asm.data_ptr<int64_t>()[k] = (int64_t)i; ctg_local.data_ptr<int64_t>()[k] = (int64_t)r;
                ctg_boff.data_ptr<int64_t>()[k] = (int64_t)X.rec[r].boff; ctg_len.data_ptr<int64_t>()[k] = (int64_t)X.rec[r].len; }
            X.attach(tp + i, compact ? nullptr : offp + bb[i], compact ? nullptr : stp + bb[i], Pp + pb[i]);
            if (compact) { X.ct = &ct; X.gb0 = (uint64_t)blk0.data_ptr<int64_t>()[i]; X.g0 = (uint64_t)grp_base.data_ptr<int64_t>()[i]; } });
    }
    // two-level block table from the archives' full tables (still owned by them at this point)
    void build_compact() {
        const int64_t n = (int64_t)A.size(); int64_t nb = 0, ng = 0;
        std::vector<int64_t> b0(n), g0(n); for (int64_t i = 0; i < n; i++) { b0[i] = nb; g0[i] = ng; nb += (int64_t)A[i]->nb; ng += (int64_t)((A[i]->nb + PV_G - 1) / PV_G); }
        grp_off = torch::empty({std::max<int64_t>(1, ng)}, torch::kInt64); grp_st = torch::empty({std::max<int64_t>(1, ng)}, torch::kInt64);
        blen = torch::empty({std::max<int64_t>(1, nb)}, torch::kInt16); bcode = torch::empty({std::max<int64_t>(1, nb)}, torch::kInt16);
        grp_base = torch::tensor(g0, torch::kInt64); blk0 = torch::tensor(b0, torch::kInt64);
        int64_t* go = grp_off.data_ptr<int64_t>(); int64_t* gs = grp_st.data_ptr<int64_t>(); uint16_t* L = (uint16_t*)blen.data_ptr<int16_t>(); int16_t* C = bcode.data_ptr<int16_t>();
        std::vector<int64_t> eb, est;
        for (int64_t i = 0; i < n; i++) { const Archive& X = *A[i];
            for (uint64_t b = 0; b < X.nb; b++) { const uint64_t len = X.off_own[b + 1] - X.off_own[b];
                if (len > 65535) throw Err(X.path + ": block payload over 65535 B - open without compact_blocks");
                L[b0[i] + b] = (uint16_t)len;
                if (b % PV_G == 0) { go[g0[i] + b / PV_G] = X.off_own[b]; gs[g0[i] + b / PV_G] = X.st_own[b]; C[b0[i] + b] = 0; continue; }
                const int64_t pc = X.st_own[b - 1] >> 1, pd = X.st_own[b - 1] & 1, pred = pd ? pc - (int64_t)Q : pc + (int64_t)Q;
                const int64_t d = (X.st_own[b] >> 1) - pred, dir = X.st_own[b] & 1;
                if (d >= -16383 && d <= 16383) C[b0[i] + b] = (int16_t)(d * 2 + dir);
                else { C[b0[i] + b] = PV_EXC; eb.push_back(b0[i] + (int64_t)b); est.push_back(X.st_own[b]); } } }
        exc_b = torch::tensor(eb.empty() ? std::vector<int64_t>{0} : eb, torch::kInt64); exc_st = torch::tensor(est.empty() ? std::vector<int64_t>{0} : est, torch::kInt64);
        ct.grp_off = go; ct.grp_st = gs; ct.len = L; ct.code = C; ct.exc_b = exc_b.data_ptr<int64_t>(); ct.exc_st = exc_st.data_ptr<int64_t>(); ct.nexc = (int64_t)eb.size(); ct.Q = Q;
    }
    const Archive& at(int64_t i) const { if (i < 0 || i >= (int64_t)A.size()) throw Err("assembly index out of range"); return *A[i]; }

    // windows: asm [n] (int), start [n] (int64, offset in the assembly's base stream), W; rc [n] bool/uint8 or None;
    // tokens: A/a 0, C/c 1, G/g 2, T/t 3, else 4 (case then irrelevant); apply_case: lower-case runs as in the FASTA
    torch::Tensor windows(torch::Tensor wasm, torch::Tensor wstart, int64_t W, c10::optional<torch::Tensor> rc, bool tokens, bool apply_case, bool verify, int64_t T) const {
        wasm = wasm.to(torch::kCPU, torch::kInt64).contiguous(); wstart = wstart.to(torch::kCPU, torch::kInt64).contiguous();
        const int64_t n = wasm.numel(); if (wstart.numel() != n) throw Err("assembly and start counts differ"); if (W < 0) throw Err("negative W");
        torch::Tensor r; if (rc && rc->defined()) { r = rc->to(torch::kCPU, torch::kUInt8).contiguous(); if (r.numel() != n) throw Err("reverse_complement mask size"); }
        const int64_t* ap = wasm.data_ptr<int64_t>(); const int64_t* sp = wstart.data_ptr<int64_t>(); const uint8_t* rp = r.defined() ? r.data_ptr<uint8_t>() : nullptr;
        for (int64_t w = 0; w < n; w++) { const Archive& X = at(ap[w]); if (sp[w] < 0 || (uint64_t)sp[w] + (uint64_t)W > X.nbases) throw Err("window outside the assembly"); }
        torch::Tensor out = torch::empty({n, W}, torch::kUInt8); uint8_t* op = out.data_ptr<uint8_t>(); std::atomic<int> err{0};
        const uint64_t nT = (uint64_t)std::max<int64_t>(1, T > 0 ? T : threads);
        // each thread its scratch: windows are dealt in contiguous ranges
        const uint64_t chunk = std::max<uint64_t>(1, ((uint64_t)n + nT * 4 - 1) / (nT * 4));
        parallel(((uint64_t)n + chunk - 1) / chunk, (int)nT, [&](uint64_t c) {
            Scratch S(Q);
            for (uint64_t w = c * chunk; w < std::min<uint64_t>((uint64_t)n, (c + 1) * chunk) && !err; w++) {
                uint8_t* o = op + w * (uint64_t)W; const int e = with_ref([&](const auto& R) { return window_cpu(*A[ap[w]], R, ref_n, (uint64_t)sp[w], (uint64_t)W, o, S, apply_case && !tokens, verify); });
                if (e) { int z = 0; err.compare_exchange_strong(z, e); break; }
                if (rp && rp[w]) { for (int64_t i = 0, j = W - 1; i < j; i++, j--) { const uint8_t t = comp_case(o[i]); o[i] = comp_case(o[j]); o[j] = t; } if (W & 1) o[W / 2] = comp_case(o[W / 2]); }
                if (tokens) for (int64_t i = 0; i < W; i++) o[i] = token(o[i]);
            } });
        if (err) throw Err(err == 2 ? "block XXH3 mismatch (corrupt archive)" : "block decode failed (corrupt archive)");
        return out;
    }
    // full decode of assembly i -> the FASTA file bytes; verify: every block XXH3 (if present) and the FASTA XXH3
    py::bytes fasta(int64_t i, bool verify) const {
        const Archive& X = at(i); uint64_t size = 0;
        for (auto& r : X.rec) size += 2 + r.hdr.size() + r.len + (r.lw ? (r.len + r.lw - 1) / r.lw : 0);
        PyObject* o = PyBytes_FromStringAndSize(nullptr, (Py_ssize_t)size); if (!o) throw py::error_already_set();
        py::bytes res = py::reinterpret_steal<py::bytes>(o); char* dst = PyBytes_AS_STRING(o); int err = 0;
        {   py::gil_scoped_release nogil;
            std::vector<uint8_t> B(X.nbases + 64); std::atomic<int> e{0};
            const uint64_t per = 64; parallel((X.nb + per - 1) / per, threads, [&](uint64_t c) { Scratch S(Q);
                for (uint64_t b = c * per; b < std::min(X.nb, (c + 1) * per) && !e; b++) { const int k = with_ref([&](const auto& R) { return decode_block(X, R, ref_n, b, &B[b * Q], S, verify); }); if (k) { int z = 0; e.compare_exchange_strong(z, k); } } });
            err = e;
            if (!err) {
                for (size_t r = 0; r < X.low_s.size(); r++) for (uint64_t x = X.low_s[r]; x < X.low_s[r] + X.low_l[r]; x++) B[x] |= 0x20;
                char* p = dst;
                for (auto& r : X.rec) { *p++ = '>'; memcpy(p, r.hdr.data(), r.hdr.size()); p += r.hdr.size(); *p++ = '\n';
                    for (uint64_t x = 0; x < r.len; x += r.lw) { const uint64_t k = std::min<uint64_t>(r.lw, r.len - x); memcpy(p, &B[r.boff + x], k); p += k; *p++ = '\n'; } }
                if ((uint64_t)(p - dst) != size) err = 4;
                else if (verify && XXH3_64bits(dst, size) != X.fasta_xxh) err = 3;
            }
        }
        if (err) throw Err(X.path + (err == 1 ? ": block decode failed" : err == 2 ? ": block XXH3 mismatch" : err == 3 ? ": FASTA XXH3 differs from the header" : ": FASTA size"));
        return res;
    }
    py::list contigs(int64_t i) const { py::list l; for (auto& r : at(i).rec) l.append(py::make_tuple(r.name, r.hdr, r.len, r.boff, r.lw)); return l; }
    py::dict info(int64_t i) const { const Archive& X = at(i); py::dict d;
        d["path"] = X.path; d["Q"] = X.Q; d["bases"] = X.nbases; d["blocks"] = X.nb; d["file_bytes"] = X.file_bytes; d["payload_bytes"] = X.payload_bytes;
        d["block_hashes"] = !X.hashes.empty(); d["fasta_xxh3"] = x16(X.fasta_xxh); d["contigs"] = X.rec.size(); d["lower_case_runs"] = X.low_s.size(); d["reference_name"] = X.refname; return d; }
    static std::string x16(uint64_t v) { char b[17]; snprintf(b, sizeof b, "%016llx", (unsigned long long)v); return b; }
    // XXH3 of the FASTA file rebuilt from an assembly's decoded base stream (case applied; e.g. the GPU's output copied
    // to the host) and the record layout - streamed, no FASTA in memory; "%016llx" as in the manifest
    std::string fasta_xxh3(int64_t i, torch::Tensor bases) const {
        const Archive& X = at(i); bases = bases.contiguous();
        if (bases.device().type() != torch::kCPU || bases.scalar_type() != torch::kUInt8 || (uint64_t)bases.numel() != X.nbases) throw Err("bases: uint8 CPU tensor of the assembly's length");
        const uint8_t* B = bases.data_ptr<uint8_t>(); XXH3_state_t* h = XXH3_createState(); XXH3_64bits_reset(h);
        std::vector<char> buf; buf.reserve(1 << 21);
        auto flush = [&](bool all) { if (all || buf.size() >= (1u << 20)) { XXH3_64bits_update(h, buf.data(), buf.size()); buf.clear(); } };
        for (auto& r : X.rec) { buf.push_back('>'); buf.insert(buf.end(), r.hdr.begin(), r.hdr.end()); buf.push_back('\n');
            for (uint64_t x = 0; x < r.len; x += r.lw) { const uint64_t k = std::min<uint64_t>(r.lw, r.len - x); buf.insert(buf.end(), B + r.boff + x, B + r.boff + x + k); buf.push_back('\n'); flush(false); } }
        flush(true); const uint64_t v = XXH3_64bits_digest(h); XXH3_freeState(h); return x16(v);
    }
    py::dict pools() const { py::dict d;
        d["P"] = P; d["pay_base"] = pay_base; d["off"] = off; d["st"] = st; d["tab"] = tab; d["blk_base"] = blk_base; d["nbases"] = nbases;
        d["low_s"] = low_s; d["low_l"] = low_l; d["low_base"] = low_base; d["low_cnt"] = low_cnt; d["ref"] = ref;
        if (packed) { d["pw"] = pw; d["pes"] = pes; d["pel"] = pel; d["peb"] = peb; d["pne"] = pr.ne; }
        if (compact) { d["grp_off"] = grp_off; d["grp_st"] = grp_st; d["blen"] = blen; d["bcode"] = bcode; d["exc_b"] = exc_b; d["exc_st"] = exc_st; d["nexc"] = ct.nexc; d["grp_base"] = grp_base; d["blk0"] = blk0; }
        d["ctg_asm"] = ctg_asm; d["ctg_local"] = ctg_local; d["ctg_boff"] = ctg_boff; d["ctg_len"] = ctg_len; return d; }
};

}  // namespace pv

#ifdef PV_WITH_CUDA
std::vector<torch::Tensor> pv_windows_cuda(torch::Tensor P, torch::Tensor pay_base, torch::Tensor off, torch::Tensor st, torch::Tensor tab, torch::Tensor blk_base,
                                           torch::Tensor nbases, torch::Tensor low_s, torch::Tensor low_l, torch::Tensor low_base, torch::Tensor low_cnt, torch::Tensor ref,
                                           int64_t Q, torch::Tensor wasm, torch::Tensor wstart, int64_t W, c10::optional<torch::Tensor> rc, bool tokens, bool apply_case);
std::vector<torch::Tensor> pv_windows_cuda2(std::vector<torch::Tensor> t, int64_t Q, int64_t ref_n, int64_t pne, int64_t nexc, bool packed, bool compact,
                                            torch::Tensor wasm, torch::Tensor wstart, int64_t W, c10::optional<torch::Tensor> rc, bool tokens, bool apply_case);
#endif

PYBIND11_MODULE(TORCH_EXTENSION_NAME, m) {
    py::register_exception<pv::Err>(m, "FormatError", PyExc_ValueError);
    py::class_<pv::Core>(m, "Core")
        .def(py::init<const std::string&, const std::vector<std::string>&, int, bool, bool>(), py::arg("reference"), py::arg("paths"), py::arg("threads") = 0,
             py::arg("packed_reference") = false, py::arg("compact_blocks") = false, py::call_guard<py::gil_scoped_release>())
        .def("windows", &pv::Core::windows, py::arg("asm"), py::arg("start"), py::arg("W"), py::arg("rc") = py::none(), py::arg("tokens") = false, py::arg("apply_case") = true,
             py::arg("verify") = false, py::arg("threads") = 0, py::call_guard<py::gil_scoped_release>())
        .def("fasta", &pv::Core::fasta, py::arg("i"), py::arg("verify") = true)
        .def("fasta_xxh3", &pv::Core::fasta_xxh3, py::call_guard<py::gil_scoped_release>())
        .def("contigs", &pv::Core::contigs).def("info", &pv::Core::info).def("pools", &pv::Core::pools)
        .def_property_readonly("Q", [](const pv::Core& c) { return c.Q; })
        .def_property_readonly("packed_reference", [](const pv::Core& c) { return c.packed; })
        .def_property_readonly("compact_blocks", [](const pv::Core& c) { return c.compact; })
        .def_property_readonly("n", [](const pv::Core& c) { return c.A.size(); })
        .def_property_readonly("ref_bases", [](const pv::Core& c) { return c.ref_n; })
        .def_property_readonly("ref_sha256", [](const pv::Core& c) { return pv::hex(c.ref_sha, 32); });
    m.def("sha256_file_bases", [](const std::string& p) { std::vector<uint8_t> b = pv::slurp(p); const uint64_t n = pv::fasta_bases_inplace(b.data(), b.size(), p); uint8_t s[32]; rr_sha256(b.data(), n, s); return pv::hex(s, 32); },
          py::call_guard<py::gil_scoped_release>());
#ifdef PV_WITH_CUDA
    m.def("windows_cuda", &pv_windows_cuda);
    m.def("windows_cuda2", &pv_windows_cuda2);
    m.attr("with_cuda") = true;
#else
    m.attr("with_cuda") = false;
#endif
}
