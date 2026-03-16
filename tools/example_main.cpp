#include "core/fields.hpp"
#include "core/grid.hpp"
#include "core/runtime.hpp"
#include "solver/momentum_terms.hpp"

#include <iostream>

int main() {
  const solver::BuildInfo build_info = solver::get_build_info();
  const solver::Grid grid{16, 8, 4, 0.1, 0.1, 0.1};
  const solver::PressureField pressure{grid};
  const solver::VelocityField velocity{grid};
  const solver::AdvectionOptions advection_options{};
  const solver::CflDiagnostics cfl = solver::compute_advective_cfl(velocity, 0.01);

  std::cout << solver::format_build_banner(build_info) << '\n';
  std::cout << "pressure_storage_cells: " << pressure.size() << '\n';
  std::cout << "u_active_nx: " << velocity.x.layout().active_extent().nx << '\n';
  std::cout << "advection: " << solver::describe(advection_options) << '\n';
  std::cout << "max_cfl: " << cfl.max_cfl << '\n';
  std::cout << "example_status: ok" << '\n';

  return build_info.supported_runtime_platform ? 0 : 1;
}
