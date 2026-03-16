#pragma once

#include "solver/momentum_terms.hpp"
#include "solver/projection.hpp"

#include <string>
#include <vector>

namespace solver {

enum class ChannelFlowCase : int {
  couette = 0,
  poiseuille = 1,
};

struct ChannelFlowConfig {
  ChannelFlowCase case_kind = ChannelFlowCase::couette;
  int nx = 128;
  int ny = 128;
  double viscosity = 0.1;
  double top_velocity = 1.0;
  double pressure_drop = 0.8;
  double cfl_limit = 0.5;
  int steps = 8;
  int poisson_max_iterations = 200;
  double poisson_tolerance = 1.0e-10;
  bool validate_profile = true;
  AdvectionOptions advection{};
};

struct ChannelFlowStepMetrics {
  int step = 0;
  double time = 0.0;
  double dt = 0.0;
  double max_cfl = 0.0;
  double max_velocity_change = 0.0;
  double divergence_l2 = 0.0;
  int pressure_iterations = 0;
  double pressure_relative_residual = 0.0;
};

struct WallNormalProfile {
  std::vector<double> coordinate;
  std::vector<double> value;
};

struct ChannelFlowValidation {
  std::string reference_dataset;
  double relative_l2_error = 0.0;
  double divergence_l2 = 0.0;
  bool pass = false;
};

struct ChannelFlowResult {
  ChannelFlowConfig config{};
  ChannelFlowStepMetrics final_step{};
  WallNormalProfile streamwise_profile{};
  ChannelFlowValidation validation{};
};

[[nodiscard]] ChannelFlowConfig default_channel_flow_config();
[[nodiscard]] ChannelFlowConfig load_channel_flow_config(const std::string& path);
[[nodiscard]] std::string to_string(ChannelFlowCase case_kind);
[[nodiscard]] std::string describe(const ChannelFlowConfig& config);
[[nodiscard]] BoundaryConditionSet make_channel_flow_boundary_conditions(
    const ChannelFlowConfig& config);
[[nodiscard]] ChannelFlowValidation validate_channel_flow_profile(const ChannelFlowResult& result);
[[nodiscard]] ChannelFlowResult run_channel_flow(const ChannelFlowConfig& config);

}  // namespace solver
