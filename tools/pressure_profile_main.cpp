#include "linsolve/poisson_solver.hpp"

#include <chrono>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

double pi() {
  return std::acos(-1.0);
}

template <typename Field, typename ValueFn>
void fill_storage(Field& field, ValueFn&& value_fn) {
  const solver::Extent3D storage = field.layout().storage_extent();

  for(int k = 0; k < storage.nz; ++k) {
    const double z = field.layout().coordinate_for_storage_index(solver::Axis::z, k);
    for(int j = 0; j < storage.ny; ++j) {
      const double y = field.layout().coordinate_for_storage_index(solver::Axis::y, j);
      for(int i = 0; i < storage.nx; ++i) {
        const double x = field.layout().coordinate_for_storage_index(solver::Axis::x, i);
        field(i, j, k) = value_fn(x, y, z);
      }
    }
  }
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const int resolution = argc > 1 ? std::stoi(argv[1]) : 256;
    const int repeats = argc > 2 ? std::stoi(argv[2]) : 24;
    if(resolution <= 0 || repeats <= 0) {
      throw std::invalid_argument("resolution and repeats must be positive");
    }

    const double spacing = 1.0 / static_cast<double>(resolution);
    const solver::Grid grid{resolution, resolution, 1, spacing, spacing, 1.0, 1};
    solver::ScalarField rhs{grid};
    solver::PressureField pressure{grid};
    solver::PressureBoundarySet boundary_conditions;
    for(auto& face : boundary_conditions.faces) {
      face.type = solver::PressureBoundaryType::dirichlet;
      face.value = 0.0;
    }

    fill_storage(rhs, [](const double x, const double y, double) {
      return 2.0 * pi() * pi() * std::sin(pi() * x) * std::sin(pi() * y);
    });

    const solver::ProjectionOptions options{
        .dt = 1.0,
        .density = 1.0,
        .poisson_max_iterations = 4000,
        .poisson_tolerance = 1.0e-10,
    };

    int iteration_sum = 0;
    double residual_sum = 0.0;
    int multigrid_levels = 0;

    const auto start = std::chrono::steady_clock::now();
    for(int repeat = 0; repeat < repeats; ++repeat) {
      pressure.fill(0.0);
      const solver::PoissonSolveDiagnostics diagnostics =
          solver::linsolve::solve_pressure_poisson(rhs, boundary_conditions, options, pressure);
      if(!diagnostics.converged) {
        throw std::runtime_error("pressure profile benchmark encountered a non-converged solve");
      }
      iteration_sum += diagnostics.iterations;
      residual_sum += diagnostics.relative_residual;
      multigrid_levels = diagnostics.multigrid_levels;
    }
    const auto stop = std::chrono::steady_clock::now();

    const double elapsed_seconds =
        std::chrono::duration_cast<std::chrono::duration<double>>(stop - start).count();
    const double average_seconds = elapsed_seconds / static_cast<double>(repeats);
    const double unknowns_per_second =
        static_cast<double>(resolution) * static_cast<double>(resolution) *
        static_cast<double>(repeats) / elapsed_seconds;

    std::cout << "kernel: pressure_poisson" << '\n';
    std::cout << "nx: " << resolution << '\n';
    std::cout << "ny: " << resolution << '\n';
    std::cout << "repeats: " << repeats << '\n';
    std::cout << "elapsed_seconds: " << elapsed_seconds << '\n';
    std::cout << "average_seconds_per_repeat: " << average_seconds << '\n';
    std::cout << "logical_cells: " << static_cast<std::size_t>(resolution) * static_cast<std::size_t>(resolution)
              << '\n';
    std::cout << "unknown_updates_per_second: " << unknowns_per_second << '\n';
    std::cout << "average_iterations: "
              << static_cast<double>(iteration_sum) / static_cast<double>(repeats) << '\n';
    std::cout << "average_relative_residual: "
              << residual_sum / static_cast<double>(repeats) << '\n';
    std::cout << "multigrid_levels: " << multigrid_levels << '\n';
    return 0;
  } catch(const std::exception& exception) {
    std::cerr << "solver_pressure_profile failed: " << exception.what() << '\n';
    return 1;
  }
}
