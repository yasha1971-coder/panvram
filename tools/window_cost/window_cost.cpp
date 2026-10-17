// SPDX-License-Identifier: MIT
// window_cost.cpp - window-cost verdict for refrel3 v1 archives (model B of the hw-apex window-law job, protocol v1.1),
// built on the public reader of this repository only (csrc/pv_v1.h: open_v1 with every check, decode_block).
//
//   window_cost truth <reference.fa> <archive dir> <names.txt> <requests.tsv> [threads]
//       every archive of names.txt decoded in full (block hashes checked when present, bases XXH3 == header), then
//       SAMPLE <name> <bases> <XXH3 hex> and REQ <W> <id> <SHA-256 of the request's upper-case bases>; stdout = truth table
//   window_cost run <reference.fa> <archive dir> <names.txt> <requests.tsv> <truth.tsv> [--verify-blocks]
//       one thread, in process. D_Q = canonical bases of the cohort / median of 9 full sequential decodes after
//       3 warm-ups (output buffer allocated outside the timer; XXH3 of every decoded assembly checked after its
//       interval). Windows: per request one timed fetch (contig lookup by name + the touched blocks decoded + copy),
//       SHA-256 checked after the timer; p50 nearest rank; c0 = p50(W = 1) - Q / D_Q; predicted
//       t(W) = c0 + (W + Q - 1) / D_Q; error = 100 (p50 - pred) / pred; PASS if c0 >= 0 and |error| <= 20 % at
//       W = 1024, 8192, 65536 and no request failed. Block hashes are not checked unless --verify-blocks.
//       WINDOW_COST_RAW=<file>: every timed request as "W id nanoseconds". Exit 0 PASS, 5 FAIL, 4 decode failure,
//       2 refused input.
// requests.tsv: W id assembly contig start0 end0 (half-open, canonical domain = upper-case bases of the contig record),
// from tools/window_cost/window_requests.py (the frozen generator of the hw-apex job, seed 20261003, 10 000 per W).
#include "pv_v1.h"
#include <atomic>
#include <chrono>
#include <cmath>
#include <map>
#include <thread>

using namespace pv;

struct Req { int W; std::string id, asmb, ctg; uint64_t s, e; };

static std::vector<std::vector<std::string>> read_tsv(const std::string& p) {
    std::vector<std::vector<std::string>> rows; FILE* f = fopen(p.c_str(), "r"); if (!f) throw Err(p + ": cannot open");
    char* line = nullptr; size_t cap = 0; ssize_t n;
    while ((n = getline(&line, &cap, f)) >= 0) { std::string s(line, (size_t)n); while (!s.empty() && (s.back() == '\n' || s.back() == '\r')) s.pop_back();
        std::vector<std::string> r; size_t a = 0; for (;;) { const size_t t = s.find('\t', a); r.push_back(s.substr(a, t == std::string::npos ? std::string::npos : t - a)); if (t == std::string::npos) break; a = t + 1; }
        rows.push_back(r); }
    free(line); fclose(f); return rows;
}
static std::vector<std::string> read_names(const std::string& p) { std::vector<std::string> v; for (auto& r : read_tsv(p)) if (!r.empty() && !r[0].empty()) v.push_back(r[0]); return v; }
static std::vector<Req> read_reqs(const std::string& p) { std::vector<Req> v; for (auto& r : read_tsv(p)) if (r.size() >= 6 && r[0] != "W") v.push_back({atoi(r[0].c_str()), r[1], r[2], r[3], strtoull(r[4].c_str(), 0, 10), strtoull(r[5].c_str(), 0, 10)}); return v; }
static std::string archive_path(const std::string& dir, const std::string& nm) { return dir + "/" + nm + ".q16k.rr3"; }
static std::string sha_hex(const uint8_t* d, size_t n) { uint8_t md[32]; rr_sha256(d, n, md); return hex(md, 32); }
static inline uint64_t now_ns() { return (uint64_t)std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now().time_since_epoch()).count(); }

struct Reference { std::vector<uint8_t> R; uint8_t sha[32]; uint64_t n = 0; };
static Reference load_reference(const std::string& path) {
    Reference r; r.R = slurp(path); r.n = fasta_bases_inplace(r.R.data(), r.R.size(), path); r.R.resize(r.n); r.R.shrink_to_fit(); rr_sha256(r.R.data(), r.n, r.sha); return r;
}
static std::unique_ptr<Archive> open_archive(const std::string& dir, const std::string& nm, const Reference& ref) {
    const std::string p = archive_path(dir, nm); const std::vector<uint8_t> f = slurp(p);
    std::unique_ptr<Archive> X(new Archive()); open_v1(f, p, ref.sha, ref.n, *X); return X;   // the payload is copied into X; f is released here
}

static int cmd_truth(int argc, char** argv) {
    const Reference ref = load_reference(argv[2]); const std::string dir = argv[3];
    const std::vector<std::string> names = read_names(argv[4]); const std::vector<Req> q = read_reqs(argv[5]); const int T = argc > 6 ? std::max(1, atoi(argv[6])) : 16;
    for (auto& nm : names) {
        std::unique_ptr<Archive> X = open_archive(dir, nm, ref); const Archive& A_ = *X;
        std::vector<uint8_t> A(A_.nbases + 64); std::atomic<uint64_t> next{0}; std::atomic<int> err{0}; std::vector<std::thread> th;
        for (int t = 0; t < T; t++) th.emplace_back([&] { Scratch S(A_.Q); for (;;) { const uint64_t b = next.fetch_add(1); if (b >= A_.nb || err.load()) break;
            if (decode_block(A_, ref.R.data(), ref.n, b, &A[b * A_.Q], S, true)) err.store(1); } });
        for (auto& x : th) x.join();
        if (err.load()) { fprintf(stderr, "%s: block decode or block hash failed\n", nm.c_str()); return 3; }
        const uint64_t xx = XXH3_64bits(A.data(), A_.nbases);
        if (xx != A_.bases_xxh) { fprintf(stderr, "%s: bases XXH3 differs from the header\n", nm.c_str()); return 3; }
        printf("SAMPLE\t%s\t%llu\t%016llx\n", nm.c_str(), (unsigned long long)A_.nbases, (unsigned long long)xx);
        std::map<std::string, uint64_t> bo; for (auto& r : A_.rec) bo[r.name] = r.boff;
        for (auto& r : q) if (r.asmb == nm) { auto it = bo.find(r.ctg); if (it == bo.end()) { fprintf(stderr, "%s: contig %s not in the archive\n", nm.c_str(), r.ctg.c_str()); return 2; }
            printf("REQ\t%d\t%s\t%s\n", r.W, r.id.c_str(), sha_hex(&A[it->second + r.s], r.e - r.s).c_str()); }
        fflush(stdout);
    }
    return 0;
}

static int cmd_run(int argc, char** argv) {
    const std::string dir = argv[3]; const std::vector<std::string> names = read_names(argv[4]); const std::vector<Req> q = read_reqs(argv[5]);
    const bool verify = argc > 7 && !strcmp(argv[7], "--verify-blocks");
    std::map<std::string, std::pair<uint64_t, uint64_t>> ts; std::map<std::pair<int, std::string>, std::string> tr;
    for (auto& r : read_tsv(argv[6])) { if (r.size() >= 4 && r[0] == "SAMPLE") ts[r[1]] = {strtoull(r[2].c_str(), 0, 10), strtoull(r[3].c_str(), 0, 16)}; else if (r.size() >= 4 && r[0] == "REQ") tr[{atoi(r[1].c_str()), r[2]}] = r[3]; }
    const uint64_t t0 = now_ns(); const Reference ref = load_reference(argv[2]); const uint8_t* G = ref.R.data(); const uint64_t Rn = ref.n;
    const size_t N = names.size(); std::vector<std::unique_ptr<Archive>> X(N); std::map<std::string, size_t> idx; uint64_t total = 0, maxb = 0; const uint32_t Q = 16384;
    for (size_t i = 0; i < N; i++) { const std::string& nm = names[i]; idx[nm] = i;
        try { X[i] = open_archive(dir, nm, ref); } catch (const Err& e) { fprintf(stderr, "%s\n", e.what()); return 2; }
        if (X[i]->Q != Q) { fprintf(stderr, "%s: Q differs\n", nm.c_str()); return 2; }
        total += X[i]->nbases; maxb = std::max(maxb, X[i]->nbases);
        if (!ts.count(nm) || ts[nm].first != X[i]->nbases) { fprintf(stderr, "%s: no truth / size differs\n", nm.c_str()); return 2; } }
    Scratch S(Q);
    auto block = [&](size_t i, uint64_t b, uint8_t* out) -> int { return decode_block(*X[i], G, Rn, b, out, S, verify); };
    printf("WINDOW_COST\tvariant v1\tassemblies %zu\tcanonical bases %llu\topen %.1f s\n", N, (unsigned long long)total, (now_ns() - t0) / 1e9); fflush(stdout);
    // D_Q
    std::vector<uint8_t> A(maxb + 64); std::vector<uint64_t> dq;
    for (int it = 0; it < 12; it++) { uint64_t sum = 0; bool ok = true;
        for (size_t i = 0; i < N; i++) { const uint64_t nb = X[i]->nb, nbs = X[i]->nbases;
            const uint64_t a = now_ns(); for (uint64_t b = 0; b < nb; b++) if (block(i, b, &A[b * Q])) { ok = false; break; } const uint64_t z = now_ns(); sum += z - a;
            if (!ok || XXH3_64bits(A.data(), nbs) != ts[names[i]].second) { printf("DQ\t%d\tFAILED %s\n", it, names[i].c_str()); return 4; } }
        printf("DQ\t%d\t%s\t%llu ns\tall %zu XXH3 ok\n", it, it < 3 ? "warmup" : "measured", (unsigned long long)sum, N); fflush(stdout); if (it >= 3) dq.push_back(sum); }
    std::vector<uint64_t> s9 = dq; std::sort(s9.begin(), s9.end()); const double med = (double)s9[4], DQ = total / (med * 1e-9);
    printf("DQ_RESULT\tmedian %.0f ns\tmin %llu\tmax %llu\tD_Q %.6f Gbases/s\n", med, (unsigned long long)s9[0], (unsigned long long)s9[8], DQ / 1e9);
    // windows
    std::vector<uint8_t> dst(65536 + 64); double c0 = 0; int fails = 0; std::string raw;
    for (int W : {1, 1024, 8192, 65536}) { std::vector<uint64_t> lat; int bad = 0;
        for (auto& r : q) { if (r.W != W) continue; const size_t i = idx.at(r.asmb); const std::vector<Contig>& rec = X[i]->rec;
            const uint64_t a = now_ns();
            size_t ri = 0; while (ri < rec.size() && rec[ri].name != r.ctg) ri++;
            if (ri == rec.size() || r.e > rec[ri].len) { bad++; continue; }
            const uint64_t lo = rec[ri].boff + r.s, hi = rec[ri].boff + r.e; size_t n = 0; bool ok = true;
            for (uint64_t b = lo / Q; b <= (hi - 1) / Q; b++) { if (block(i, b, S.blk.data())) { ok = false; break; } const uint64_t bs = b * Q, x0 = std::max(lo, bs), x1 = std::min<uint64_t>(hi, bs + Q); memcpy(&dst[n], S.blk.data() + x0 - bs, x1 - x0); n += x1 - x0; }
            const uint64_t z = now_ns();
            if (!ok || n != r.e - r.s || sha_hex(dst.data(), n) != tr[{W, r.id}]) { bad++; continue; }
            lat.push_back(z - a); raw += std::to_string(W) + "\t" + r.id + "\t" + std::to_string(z - a) + "\n"; }
        std::vector<uint64_t> s = lat; std::sort(s.begin(), s.end()); auto nr = [&](int p) { return s.empty() ? 0.0 : s[(p * s.size() + 99) / 100 - 1] / 1000.0; };
        const double p50 = nr(50); if (W == 1) c0 = p50 - Q / DQ * 1e6;
        const double pred = c0 + (W + Q - 1) / DQ * 1e6, err = 100 * (p50 - pred) / pred; const bool verdict = W != 1;
        const bool pass = bad == 0 && (!verdict || (c0 >= 0 && std::fabs(err) <= 20.0)); if (!pass) fails++;
        printf("WIN\tW %d\tverified %zu\tbad %d\tp50 %.3f us\tp95 %.3f\tp99 %.3f\tpred %.3f us\terror %+.2f %%\t%s\n", W, lat.size(), bad, p50, nr(95), nr(99), pred, err,
               W == 1 ? "c0-probe" : pass ? "PASS" : "FAIL"); fflush(stdout); }
    printf("RESULT\tvariant v1\tQ %u\tD_Q %.6f Gbases/s\tc0 %.3f us\tverdict %s\n", Q, DQ / 1e9, c0, fails ? "FAIL" : "PASS");
    if (const char* p = getenv("WINDOW_COST_RAW")) { FILE* f = fopen(p, "w"); if (f) { fwrite(raw.data(), 1, raw.size(), f); fclose(f); } }
    return fails ? 5 : 0;
}

int main(int argc, char** argv) {
    try {
        if (argc >= 6 && !strcmp(argv[1], "truth")) return cmd_truth(argc, argv);
        if (argc >= 7 && !strcmp(argv[1], "run")) return cmd_run(argc, argv);
    } catch (const std::exception& e) { fprintf(stderr, "%s\n", e.what()); return 2; }
    fprintf(stderr, "usage: window_cost truth <reference.fa> <archive dir> <names.txt> <requests.tsv> [threads]\n"
                    "       window_cost run <reference.fa> <archive dir> <names.txt> <requests.tsv> <truth.tsv> [--verify-blocks]\n");
    return 1;
}
