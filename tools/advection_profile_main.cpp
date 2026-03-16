#include "solver/momentum_terms.hpp"

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

std::size_t active_count(const solver::StructuredField& field) {
  return field.layout().active_range().extent().cell_count();
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const int resolution = argc > 1 ? std::stoi(argv[1]) : 256;
    const int repeats = argc > 2 ? std::stoi(argv[2]) : 200;
    if(resolution <= 0 || repeats <= 0) {
      throw std::invalid_argument("resolution and repeats must be positive");
    }

    const double domain_length = 2.0 * pi();
    const double spacing = domain_length / static_cast<double>(resolution);
    const solver::Grid grid{resolution, resolution, 1, spacing, spacing, 1.0, 1};

    solver::VelocityField velocity{grid};
    solver::VelocityField advection{grid};
    const solver::AdvectionOptions options{};

    fill_storage(velocity.x, [](const double x, const double y, double) {
      return -std::cos(x) * std::sin(y);
    });
    fill_storage(velocity.y, [](const double x, const double y, double) {
      return std::sin(x) * std::cos(y);
    });
    fill_storage(velocity.z, [](double, double, double) {
      return 0.0;
    });

    const std::size_t logical_cells =
        static_cast<std::size_t>(resolution) * static_cast<std::size_t>(resolution);
    const std::size_t lower_bound_bytes_per_repeat =
        8ull * (active_count(velocity.x) + active_count(velocity.y) + active_count(velocity.z) +
                active_count(advection.x) + active_count(advection.y) + active_count(advection.z));

    const auto start = std::chrono::steady_clock::now();
    for(int repeat = 0; repeat < repeats; ++repeat) {
      solver::compute_advection_term(velocity, options, advection);
    }
    const auto stop = std::chrono::steady_clock::now();

    const double elapsed_seconds =
        std::chrono::duration_cast<std::chrono::duration<double>>(stop - start).count();
    const double average_seconds = elapsed_seconds / static_cast<double>(repeats);
    const double cell_updates_per_second =
        static_cast<double>(logical_cells) * static_cast<double>(repeats) / elapsed_seconds;
    const double lower_bound_bandwidth_gbps =
        static_cast<double>(lower_bound_bytes_per_repeat) * static_cast<double>(repeats) /
        elapsed_seconds / 1.0e9;

    std::cout << "kernel: advection" << '\n';
    std::cout << "nx: " << resolution << '\n';
    std::cout << "ny: " << resolution << '\n';
    std::cout << "repeats: " << repeats << '\n';
    std::cout << "elapsed_seconds: " << elapsed_seconds << '\n';
    std::cout << "average_seconds_per_repeat: " << average_seconds << '\n';
    std::cout << "logical_cells: " << logical_cells << '\n';
    std::cout << "cell_updates_per_second: " << cell_updates_per_second << '\n';
    std::cout << "lower_bound_bandwidth_gbps: " << lower_bound_bandwidth_gbps << '\n';
    std::cout << "advection: " << solver::describe(options) << '\n';
    return 0;
  } catch(const std::exception& exception) {
    std::cerr << "solver_advection_profile failed: " << exception.what() << '\n';
    return 1;
  }
}
