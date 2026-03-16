#pragma once

#include "solver/momentum_terms.hpp"
#include "solver/projection.hpp"

#include <string>
#include <vector>

namespace solver {

struct LidDrivenCavityConfig {
  int nx = 128;
  int ny = 128;
  double reynolds = 100.0;
  double lid_velocity = 1.0;
  double cfl_limit = 0.5;
  int max_steps = 8000;
  int min_steps = 200;
  double steady_tolerance = 1.0e-7;
  int poisson_max_iterations = 200;
  double poisson_tolerance = 1.0e-10;
  bool validate_reference = true;
  AdvectionOptions advection{};
};

struct SimulationStepMetrics {
  int step = 0;
  double time = 0.0;
  double dt = 0.0;
  double max_cfl = 0.0;
  double max_velocity_change = 0.0;
  double divergence_l2 = 0.0;
  int pressure_iterations = 0;
  double pressure_relative_residual = 0.0;
};

struct CenterlineProfile {
  std::vector<double> coordinate;
  std::vector<double> value;
};

struct CenterlineExtrema {
  double u_vertical_max = 0.0;
  double u_vertical_min = 0.0;
  double v_horizontal_max = 0.0;
  double v_horizontal_min = 0.0;
};

struct LidDrivenCavityReference {
  std::string dataset;
  double u_vertical_max_y = 0.0;
  double u_vertical_max = 0.0;
  double u_vertical_min_y = 0.0;
  double u_vertical_min = 0.0;
  double v_horizontal_max_x = 0.0;
  double v_horizontal_max = 0.0;
  double v_horizontal_min_x = 0.0;
  double v_horizontal_min = 0.0;
};

struct LidDrivenCavityValidation {
  std::string reference_dataset;
  double u_vertical_max_sample = 0.0;
  double u_vertical_max_relative_error = 0.0;
  double u_vertical_min_sample = 0.0;
  double u_vertical_min_relative_error = 0.0;
  double v_horizontal_max_sample = 0.0;
  double v_horizontal_max_relative_error = 0.0;
  double v_horizontal_min_sample = 0.0;
  double v_horizontal_min_relative_error = 0.0;
  double divergence_l2 = 0.0;
  bool pass = false;
};

struct LidDrivenCavityResult {
  LidDrivenCavityConfig config{};
  SimulationStepMetrics final_step{};
  CenterlineProfile u_vertical_centerline{};
  CenterlineProfile v_horizontal_centerline{};
  CenterlineExtrema extrema{};
  LidDrivenCavityValidation validation{};
};

[[nodiscard]] LidDrivenCavityConfig default_lid_driven_cavity_config();
[[nodiscard]] LidDrivenCavityConfig load_lid_driven_cavity_config(const std::string& path);
[[nodiscard]] std::string describe(const LidDrivenCavityConfig& config);
[[nodiscard]] BoundaryConditionSet make_lid_driven_cavity_boundary_conditions(
    const LidDrivenCavityConfig& config);
[[nodiscard]] LidDrivenCavityReference ghia_re100_reference();
[[nodiscard]] LidDrivenCavityValidation validate_lid_driven_cavity_re100(
    const LidDrivenCavityResult& result);
[[nodiscard]] LidDrivenCavityResult run_lid_driven_cavity(const LidDrivenCavityConfig& config);

}  // namespace solver
