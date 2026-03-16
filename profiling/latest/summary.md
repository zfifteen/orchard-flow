# Performance Summary

## Hardware

- CPU: Apple M1 Max
- Performance cores: 8
- Efficiency cores: 2
- Memory: 32.0 GiB

## Findings

- Recommended default execution mode: benchmark build profile with the default unclamped scheduler policy. It was the fastest measured policy on this machine.
- Current compute-thread recommendation: 1. The solver path remains single-threaded, so thread scaling is a baseline-only study for now.
- Advection microbenchmark throughput: 1.270070e+07 cell updates/s with a 0.814 GB/s lower-bound bandwidth estimate.
- Pressure microbenchmark throughput: 2.951090e+05 unknown updates/s with 16.00 average iterations.
- End-to-end Taylor-Green benchmark throughput: 2.907983e+05 cells/s.
- Default-policy hotspot categories were dominated by pressure-solve, predictor/ADI, and advection work, with the top sampled category counts: {'field_layout': 151, 'pressure_solve': 51, 'other': 29, 'advection': 16}.

## Notes

- QoS clamp measurements come from default, utility, and background taskpolicy launches.
- Core-class shares come from xctrace Time Profiler samples on Apple Silicon and are reported as sampled P-core vs E-core shares.
- The thread-scaling section is intentionally explicit that the current solver path is single-threaded; Milestone 11 is where that changes.
