# Apple-Native Solver

This repository is the beginning of a very specific kind of CFD project: an incompressible flow solver that is meant to feel native to Apple Silicon rather than merely portable to it. The long-term goal is not to assemble a generic research code that happens to compile on a Mac. The goal is to build a solver that is numerically serious, reproducible, restartable, profileable, and understandable, while also being unapologetically tuned for the realities of Apple hardware and toolchains.

That combination matters here. A lot of scientific software grows by accumulating features first and worrying about rigor later. This project is trying to do the opposite. The solver is being designed from the outset around a clear numerical method, explicit validation gates, deterministic execution rules, and a roadmap that only moves forward when the previous layer is correct. If the end result works, it should be the kind of codebase a new engineer can open, reason about, benchmark, and trust.

Right now, the repository is still early. It does not yet contain the full Navier-Stokes timestepper described in the technical spec. What it does contain is the locked-down foundation for that system: the build environment, the project structure, the profile split between validation and benchmarking, the structured-grid and field-storage layer, the first discrete operators, a smoke-test executable, a focused test harness, and a tiny profiling helper. In other words, this repo is already opinionated about how the solver should be built and verified before the time-integration and linear-solver machinery arrives.

If you are reading this as a developer, the shortest useful summary is: this project is building toward a production-grade, Apple-Silicon-native incompressible flow solver, and the repository currently reflects Milestone 2 of that plan.

## Current Status

The repository is currently at **Milestone 2: Discrete Operators**.

Implemented today:

- Apple-Silicon-only CMake configuration using Apple Clang
- two labeled build profiles: `deterministic` and `benchmark`
- validated `Grid` representation for structured Cartesian domains
- MAC-aware pressure, scalar, and velocity field storage types
- flat contiguous `double` buffers with explicit ghost-cell-aware indexing
- boundary slab and ghost-layer helpers for future BC work
- second-order discrete `gradient`, `divergence`, and `laplacian` operators
- manufactured-solution convergence checks for the operator layer
- a smoke-test executable that reports build/runtime metadata
- a minimal test executable wired into CTest
- a simple time-based profiling helper script
- scaffolded module layout for later solver milestones beyond the infrastructure layer

What is not implemented yet:

- advection, diffusion, and projection steps
- Poisson / multigrid solver infrastructure
- checkpointing, benchmark cases, and validation harnesses beyond the current infrastructure tests

## Project Goals

The project is working toward a solver that can:

- solve the incompressible Navier-Stokes equations on structured Cartesian grids
- produce benchmark-correct results on standard CFD validation cases
- run efficiently on Apple M1 Max class hardware
- produce deterministic reference results under a locked validation profile
- support restartable long-running simulations
- expose diagnostics and profiling information instead of treating performance as an afterthought
- remain maintainable for engineers who did not originally author it

## Platform and Scope

This project is intentionally narrow in scope.

- Supported platform: macOS on Apple Silicon, with Apple Clang and CMake
- Primary hardware target: Apple M1 Max
- Primary execution path: CPU-first, multithreaded, Apple-native
- Initial mesh scope: structured Cartesian grids
- Initial physics scope: incompressible, constant-density flow

Explicitly out of scope for the v1 solver:

- x86, Linux, and Windows support
- CUDA, AMD GPU paths, and portability layers
- unstructured meshes
- adaptive mesh refinement
- compressible flow
- turbulence models
- multiphase flow

This narrow scope is deliberate. The project is optimizing for correctness and a clean implementation path before it optimizes for breadth.

## Numerical Direction

The full solver architecture is described in [TECH-SPEC.md](TECH-SPEC.md), but the current direction is already fixed:

- discretization: finite volume method on a staggered MAC grid
- convective form: conservative `∇·(u⊗u)` form
- default advection scheme: bounded second-order TVD with the `van Leer` limiter
- time integration: second-order semi-implicit pressure-correction scheme
- explicit term startup: Forward Euler bootstrap into Adams-Bashforth 2
- diffusion treatment: Crank-Nicolson with ADI Helmholtz line solves
- pressure solve: MGPCG with fixed V-cycle multigrid and damped-Jacobi smoothing
- precision policy: `double` for solution state and pressure-solver reductions
- validation default: advective CFL `<= 0.5` unless a benchmark case says otherwise

Parts of that numerical path are now implemented at the operator layer, and the rest remains the governing design contract for the upcoming milestones.

## Repository Layout

The repository is organized around the milestone structure defined in the roadmap:

```text
solver/
  core/
  operators/
  solver/
  linsolve/
  bc/
  io/
  tests/
  benchmarks/
  tools/
```

What those directories mean in practice right now:

- `core/`: runtime/build metadata plus the first structured-grid and field-storage layer
- `operators/`: second-order structured-grid discrete operators
- `tools/`: the smoke-test executable and profiling helper
- `tests/`: the minimal CTest-backed test executable
- `solver/`, `linsolve/`, `bc/`, `io/`, `benchmarks/`: scaffolded directories reserved for later milestones

The technical and roadmap documents live at the repository root and currently act as the primary design references.

## Build Profiles

The build is intentionally split into two named profiles:

- `deterministic`: the reference build for validation and reproducibility work
- `benchmark`: the performance-oriented build for clearly labeled experiments

At the current stage, those profiles differ mainly in floating-point behavior:

- `deterministic` disables `-ffast-math`
- `benchmark` enables `-ffast-math` when the compiler accepts it

Both profiles target Apple Silicon and require Apple Clang.

## Building

Configure and build the deterministic profile:

```bash
cmake --preset deterministic
cmake --build build/deterministic
```

Configure and build the benchmark profile:

```bash
cmake --preset benchmark
cmake --build build/benchmark
```

The current build produces:

- `solver_core`: a small core library for runtime/build metadata and mesh/field infrastructure
- `solver_operators`: the discrete-operator library built on top of the core field layer
- `solver_example`: a smoke-test executable under `build/<profile>/tools/`
- `solver_tests`: a minimal test executable under `build/<profile>/tests/`

## Testing

Run the deterministic-profile tests:

```bash
ctest --test-dir build/deterministic --output-on-failure
```

Run the benchmark-profile tests:

```bash
ctest --test-dir build/benchmark --output-on-failure
```

Today’s test coverage is still intentionally focused, but it now covers the full Milestone 2 gate. The current test executable checks that:

- the build profile is one of the locked profile names
- the runtime platform is Apple Silicon
- the generated build banner includes the active profile
- grid dimensions, spacings, and coordinate helpers behave as expected
- pressure-field indexing matches the flat-buffer mapping contract
- ghost layers and boundary slab helpers address the expected storage regions
- storage is aligned and unit-stride in the `i` direction
- MAC-grid cell-center and face-center placement is correct
- field storage uses `double`
- manufactured-solution error norms for gradient, divergence, and Laplacian decrease under refinement
- observed operator convergence is at least second-order in the Milestone 2 validation case

That is enough for Milestone 2. It is not meant to stand in for the full benchmark and regression suite described in the technical spec.

## Profiling

The repository currently includes a lightweight time-based profiling helper:

```bash
tools/profile_example.sh build/deterministic
```

That script runs the smoke-test executable through `/usr/bin/time -lp`, which gives a small but useful baseline for the locked environment. The roadmap also assumes future profiling with Instruments once real kernels exist.

## Roadmap

The implementation plan is spelled out in [EXECUTION_ROADMAP_V1.md](EXECUTION_ROADMAP_V1.md). The short version looks like this:

1. Milestone 0: lock the environment, build system, and validation scaffolding
2. Milestone 1: add grid, field, indexing, and ghost-cell infrastructure
3. Milestone 2: implement core discrete operators
4. Milestone 3: implement advection and diffusion terms
5. Milestone 4: implement the projection step
6. Milestone 5: implement the pressure linear solver system
7. Milestone 6 and beyond: benchmark validation, boundary-condition generalization, restart/output, verification, profiling, optimization, 3D support, and conditional Metal acceleration

The repository has completed the first three items in that sequence and is set up to move into advection and diffusion work next.

The important project rule is simple: **do not advance to the next milestone unless the current validation gate passes**.

## Reference Documents

The two main design documents in this repository are:

- [TECH-SPEC.md](TECH-SPEC.md): the numerical, architectural, determinism, and validation contract for the solver
- [EXECUTION_ROADMAP_V1.md](EXECUTION_ROADMAP_V1.md): the milestone-by-milestone implementation sequence and validation gates

If you need to understand why the repository is structured the way it is, start with those two files.

## License

Licensing has not been finalized yet. Until a license file is added to the repository, do not assume one.
