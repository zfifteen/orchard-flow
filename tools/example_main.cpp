#include "core/fields.hpp"
#include "core/grid.hpp"
#include "core/runtime.hpp"
#include "solver/momentum_terms.hpp"
#include "solver/projection.hpp"

#include <iostream>

int main() {
  const solver::BuildInfo build_info = solver::get_build_info();
  const solver::Grid grid{16, 8, 1, 0.1, 0.1, 1.0};
  const solver::PressureField pressure{grid};
  const solver::VelocityField velocity{grid};
  const solver::AdvectionOptions advection_options{};
  const solver::CflDiagnostics cfl = solver::compute_advective_cfl(velocity, 0.01);
  const solver::BoundaryConditionSet boundary_conditions = solver::BoundaryConditionSet::cavity();
  const solver::ProjectionOptions projection_options{
      .dt = 0.01,
      .density = 1.0,
      .poisson_max_iterations = 2000,
      .poisson_tolerance = 1.0e-12,
  };
  solver::PressureField projection_pressure{grid};
  solver::VelocityField corrected_velocity{grid};
  const solver::ProjectionDiagnostics projection = solver::project_velocity(
      velocity, boundary_conditions, projection_options, projection_pressure, corrected_velocity);

  std::cout << solver::format_build_banner(build_info) << '\n';
  std::cout << "pressure_storage_cells: " << pressure.size() << '\n';
  std::cout << "u_active_nx: " << velocity.x.layout().active_extent().nx << '\n';
  std::cout << "advection: " << solver::describe(advection_options) << '\n';
  std::cout << "max_cfl: " << cfl.max_cfl << '\n';
  std::cout << "projection_iterations: " << projection.pressure_solve.iterations << '\n';
  std::cout << "projection_divergence_l2_after: " << projection.divergence_l2_after << '\n';
  std::cout << "example_status: ok" << '\n';

  return build_info.supported_runtime_platform ? 0 : 1;
}
