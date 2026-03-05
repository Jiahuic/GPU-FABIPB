# Repository Guidelines

## Project Structure & Module Organization
This project is transitioning from a flat layout to a standard structure to support GPU work.
- `src/`: core solver code (`fmm`, `kernel`, `solver`, `io`, `main`).
- `include/`: public/shared headers (types, globals, kernel interfaces).
- `tests/`: regression drivers and reference checks.
- `data/`: sample inputs (move current `test_proteins/` here over time).
- `scripts/`: profiling, benchmark, and transfer helpers.
- `docs/`: design notes (GPU milestones, benchmark protocol).

Until migration is complete, root-level legacy files are valid; place new modules in target directories and keep compatibility wrappers small.

## Build, Test, and Development Commands
- `make coulomb`: current baseline build (CPU).
- `make clean && make coulomb`: clean rebuild for warning checks.
- `./coulomb <input>`: run solver locally (CPU validation).
- `make profile` (planned): stage timing report (P2P, P2M, M2M, M2L, L2L, L2P).
- `make gpu-stub` (planned): compile GPU interface stubs on non-GPU machines.

Keep CPU path always buildable so development can continue without accelerator hardware.

## Coding Style & Naming Conventions
- Language: C for baseline; CUDA/HIP files should use explicit backend suffixes (`*_cuda.cu`, `*_hip.cpp`).
- Keep legacy naming for touched code; use descriptive module names for new files (`m2l_gpu.cu`, `p2p_gpu.cu`).
- Prefer typed function pointers/prototypes in headers; avoid old-style `void (*)()` declarations.
- Minimize global state expansion; pass context structs for new GPU paths.

## Testing Guidelines
- No unit-test framework yet; use deterministic regression runs.
- Validate CPU results first, then compare GPU-machine runs against CPU references.
- Track: residual/error metrics, runtime breakdown, and solver iteration counts.
- For GPU PRs, include at least one small and one medium biomolecular case.

## Commit & Pull Request Guidelines
- Commit subjects: short imperative (`refactor: move kernel code to src/kernel`).
- Keep commits atomic by stage: `layout`, `api`, `gpu-stub`, `kernel-port`, `benchmark`.
- PRs must include: commands run, accuracy deltas, timing deltas, and hardware used.

## Cross-Machine Workflow (No Local GPU)
- Develop and verify CPU + GPU stubs locally.
- Commit all changes before transfer; avoid untracked benchmark artifacts.
- Push branch to remote, pull on GPU machine, run GPU benchmarks, then push results back.
