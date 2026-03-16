#include "solver/lid_driven_cavity.hpp"

#include "operators/discrete_operators.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace solver {

namespace {

std::string trim(const std::string& value) {
  const std::size_t first = value.find_first_not_of(" \t\r\n");
  if(first == std::string::npos) {
    return "";
  }

  const std::size_t last = value.find_last_not_of(" \t\r\n");
  return value.substr(first, last - first + 1);
}

bool parse_bool(const std::string& value) {
  if(value == "true" || value == "1" || value == "yes" || value == "on") {
    return true;
  }
  if(value == "false" || value == "0" || value == "no" || value == "off") {
    return false;
  }
  throw std::runtime_error("invalid boolean value: " + value);
}

void axpy_active(StructuredField& destination, const StructuredField& source, const double scale) {
  const IndexRange3D active = destination.layout().active_range();

  for(int k = active.k_begin; k < active.k_end; ++k) {
    for(int j = active.j_begin; j < active.j_end; ++j) {
      for(int i = active.i_begin; i < active.i_end; ++i) {
        destination(i, j, k) += scale * source(i, j, k);
      }
    }
  }
}

void axpy_velocity(VelocityField& destination, const VelocityField& source, const double scale) {
  axpy_active(destination.x, source.x, scale);
  axpy_active(destination.y, source.y, scale);
  axpy_active(destination.z, source.z, scale);
}

double max_velocity_change(const VelocityField& left, const VelocityField& right) {
  double change = 0.0;

  const auto accumulate = [&change](const StructuredField& lhs, const StructuredField& rhs) {
    const IndexRange3D active = lhs.layout().active_range();
    for(int k = active.k_begin; k < active.k_end; ++k) {
      for(int j = active.j_begin; j < active.j_end; ++j) {
        for(int i = active.i_begin; i < active.i_end; ++i) {
          change = std::max(change, std::abs(lhs(i, j, k) - rhs(i, j, k)));
        }
      }
    }
  };

  accumulate(left.x, right.x);
  accumulate(left.y, right.y);
  accumulate(left.z, right.z);
  return change;
}

void validate_config(const LidDrivenCavityConfig& config) {
  if(config.nx <= 0 || config.ny <= 0) {
    throw std::invalid_argument("lid-driven cavity config expects positive nx and ny");
  }
  if(config.reynolds <= 0.0) {
    throw std::invalid_argument("lid-driven cavity config expects reynolds > 0");
  }
  if(config.lid_velocity <= 0.0) {
    throw std::invalid_argument("lid-driven cavity config expects lid_velocity > 0");
  }
  if(config.cfl_limit <= 0.0) {
    throw std::invalid_argument("lid-driven cavity config expects cfl_limit > 0");
  }
  if(config.max_steps <= 0 || config.min_steps < 0) {
    throw std::invalid_argument("lid-driven cavity config expects non-negative step controls");
  }
  if(config.poisson_max_iterations <= 0 || config.poisson_tolerance <= 0.0) {
    throw std::invalid_argument("lid-driven cavity config expects positive Poisson controls");
  }
}

double fixed_dt(const LidDrivenCavityConfig& config, const Grid& grid) {
  const double denominator = config.lid_velocity / grid.dx + config.lid_velocity / grid.dy;
  return config.cfl_limit / denominator;
}

CenterlineProfile sample_vertical_centerline_u(const VelocityField& velocity) {
  const FieldLayout& layout = velocity.x.layout();
  const IndexRange3D active = layout.active_range();
  const int i_center = active.i_begin + layout.active_extent().nx / 2;

  CenterlineProfile profile;
  profile.coordinate.reserve(static_cast<std::size_t>(layout.active_extent().ny));
  profile.value.reserve(static_cast<std::size_t>(layout.active_extent().ny));

  for(int j = active.j_begin; j < active.j_end; ++j) {
    profile.coordinate.push_back(layout.coordinate_for_storage_index(Axis::y, j));
    profile.value.push_back(velocity.x(i_center, j, active.k_begin));
  }

  return profile;
}

CenterlineProfile sample_horizontal_centerline_v(const VelocityField& velocity) {
  const FieldLayout& layout = velocity.y.layout();
  const IndexRange3D active = layout.active_range();
  const int j_center = active.j_begin + layout.active_extent().ny / 2;

  CenterlineProfile profile;
  profile.coordinate.reserve(static_cast<std::size_t>(layout.active_extent().nx));
  profile.value.reserve(static_cast<std::size_t>(layout.active_extent().nx));

  for(int i = active.i_begin; i < active.i_end; ++i) {
    profile.coordinate.push_back(layout.coordinate_for_storage_index(Axis::x, i));
    profile.value.push_back(velocity.y(i, j_center, active.k_begin));
  }

  return profile;
}

CenterlineExtrema compute_centerline_extrema(const CenterlineProfile& u_profile,
                                             const CenterlineProfile& v_profile) {
  if(u_profile.value.size() < 3 || v_profile.value.size() < 3) {
    throw std::invalid_argument("centerline extrema require at least three samples per profile");
  }

  const auto u_begin = u_profile.value.begin() + 1;
  const auto u_end = u_profile.value.end() - 1;
  const auto v_begin = v_profile.value.begin() + 1;
  const auto v_end = v_profile.value.end() - 1;
  const auto [u_min_it, u_max_it] = std::minmax_element(u_begin, u_end);
  const auto [v_min_it, v_max_it] = std::minmax_element(v_begin, v_end);

  return CenterlineExtrema{
      .u_vertical_max = *u_max_it,
      .u_vertical_min = *u_min_it,
      .v_horizontal_max = *v_max_it,
      .v_horizontal_min = *v_min_it,
  };
}

double relative_error(const double value, const double reference) {
  return std::abs(value - reference) / std::abs(reference);
}

double interpolate_profile_value(const CenterlineProfile& profile, const double coordinate) {
  if(profile.coordinate.size() != profile.value.size() || profile.coordinate.size() < 2) {
    throw std::invalid_argument("centerline profile interpolation requires matching sample arrays");
  }
  if(coordinate < profile.coordinate.front() || coordinate > profile.coordinate.back()) {
    throw std::out_of_range("reference coordinate falls outside the sampled centerline range");
  }

  const auto upper = std::lower_bound(profile.coordinate.begin(), profile.coordinate.end(), coordinate);
  if(upper == profile.coordinate.begin()) {
    return profile.value.front();
  }
  if(upper == profile.coordinate.end()) {
    return profile.value.back();
  }
  if(*upper == coordinate) {
    return profile.value[static_cast<std::size_t>(upper - profile.coordinate.begin())];
  }

  const std::size_t upper_index = static_cast<std::size_t>(upper - profile.coordinate.begin());
  const std::size_t lower_index = upper_index - 1;
  if(profile.coordinate.size() >= 4) {
    const std::size_t start =
        std::min(lower_index > 0 ? lower_index - 1 : 0, profile.coordinate.size() - 4);
    double interpolated = 0.0;
    for(std::size_t i = 0; i < 4; ++i) {
      const std::size_t node = start + i;
      double basis = 1.0;
      for(std::size_t j = 0; j < 4; ++j) {
        if(i == j) {
          continue;
        }
        const std::size_t other = start + j;
        basis *= (coordinate - profile.coordinate[other]) /
                 (profile.coordinate[node] - profile.coordinate[other]);
      }
      interpolated += basis * profile.value[node];
    }
    return interpolated;
  }

  const double lower_coordinate = profile.coordinate[lower_index];
  const double upper_coordinate = profile.coordinate[upper_index];
  const double weight = (coordinate - lower_coordinate) / (upper_coordinate - lower_coordinate);
  return (1.0 - weight) * profile.value[lower_index] + weight * profile.value[upper_index];
}

}  // namespace

LidDrivenCavityConfig default_lid_driven_cavity_config() {
  return LidDrivenCavityConfig{};
}

LidDrivenCavityConfig load_lid_driven_cavity_config(const std::string& path) {
  std::ifstream input(path);
  if(!input.is_open()) {
    throw std::runtime_error("unable to open cavity config: " + path);
  }

  LidDrivenCavityConfig config = default_lid_driven_cavity_config();
  std::string line;
  int line_number = 0;
  while(std::getline(input, line)) {
    ++line_number;
    const std::size_t comment = line.find('#');
    if(comment != std::string::npos) {
      line = line.substr(0, comment);
    }

    line = trim(line);
    if(line.empty()) {
      continue;
    }

    const std::size_t separator = line.find('=');
    if(separator == std::string::npos) {
      throw std::runtime_error("invalid config line " + std::to_string(line_number) + ": " + line);
    }

    const std::string key = trim(line.substr(0, separator));
    const std::string value = trim(line.substr(separator + 1));

    if(key == "nx") {
      config.nx = std::stoi(value);
    } else if(key == "ny") {
      config.ny = std::stoi(value);
    } else if(key == "reynolds") {
      config.reynolds = std::stod(value);
    } else if(key == "lid_velocity") {
      config.lid_velocity = std::stod(value);
    } else if(key == "cfl_limit") {
      config.cfl_limit = std::stod(value);
    } else if(key == "max_steps") {
      config.max_steps = std::stoi(value);
    } else if(key == "min_steps") {
      config.min_steps = std::stoi(value);
    } else if(key == "steady_tolerance") {
      config.steady_tolerance = std::stod(value);
    } else if(key == "poisson_max_iterations") {
      config.poisson_max_iterations = std::stoi(value);
    } else if(key == "poisson_tolerance") {
      config.poisson_tolerance = std::stod(value);
    } else if(key == "validate_reference") {
      config.validate_reference = parse_bool(value);
    } else {
      throw std::runtime_error("unsupported cavity config key: " + key);
    }
  }

  validate_config(config);
  return config;
}

std::string describe(const LidDrivenCavityConfig& config) {
  std::ostringstream builder;
  builder << "nx=" << config.nx << ", ny=" << config.ny << ", Re=" << config.reynolds
          << ", lid_velocity=" << config.lid_velocity << ", cfl_limit=" << config.cfl_limit
          << ", max_steps=" << config.max_steps << ", min_steps=" << config.min_steps
          << ", steady_tolerance=" << config.steady_tolerance
          << ", poisson_max_iterations=" << config.poisson_max_iterations
          << ", poisson_tolerance=" << config.poisson_tolerance
          << ", validate_reference=" << (config.validate_reference ? "true" : "false")
          << ", advection=" << describe(config.advection);
  return builder.str();
}

BoundaryConditionSet make_lid_driven_cavity_boundary_conditions(
    const LidDrivenCavityConfig& config) {
  BoundaryConditionSet boundary_conditions = BoundaryConditionSet::cavity();
  boundary_conditions[BoundaryFace::y_max].type = PhysicalBoundaryType::prescribed_velocity;
  boundary_conditions[BoundaryFace::y_max].velocity = {config.lid_velocity, 0.0, 0.0};
  boundary_conditions[BoundaryFace::z_min].type = PhysicalBoundaryType::symmetry;
  boundary_conditions[BoundaryFace::z_max].type = PhysicalBoundaryType::symmetry;
  return boundary_conditions;
}

LidDrivenCavityReference ghia_re100_reference() {
  return LidDrivenCavityReference{
      .dataset = "Ghia1982_Re100_centerline_extrema",
      .u_vertical_max_y = 0.9766,
      .u_vertical_max = 0.84123,
      .u_vertical_min_y = 0.4531,
      .u_vertical_min = -0.21090,
      .v_horizontal_max_x = 0.2344,
      .v_horizontal_max = 0.17527,
      .v_horizontal_min_x = 0.8047,
      .v_horizontal_min = -0.24533,
  };
}

LidDrivenCavityValidation validate_lid_driven_cavity_re100(const LidDrivenCavityResult& result) {
  const LidDrivenCavityReference reference = ghia_re100_reference();
  const double u_vertical_max_sample =
      result.u_vertical_centerline.coordinate.empty()
          ? result.extrema.u_vertical_max
          : interpolate_profile_value(result.u_vertical_centerline, reference.u_vertical_max_y);
  const double u_vertical_min_sample =
      result.u_vertical_centerline.coordinate.empty()
          ? result.extrema.u_vertical_min
          : interpolate_profile_value(result.u_vertical_centerline, reference.u_vertical_min_y);
  const double v_horizontal_max_sample =
      result.v_horizontal_centerline.coordinate.empty()
          ? result.extrema.v_horizontal_max
          : interpolate_profile_value(result.v_horizontal_centerline, reference.v_horizontal_max_x);
  const double v_horizontal_min_sample =
      result.v_horizontal_centerline.coordinate.empty()
          ? result.extrema.v_horizontal_min
          : interpolate_profile_value(result.v_horizontal_centerline, reference.v_horizontal_min_x);
  LidDrivenCavityValidation validation{
      .reference_dataset = reference.dataset,
      .u_vertical_max_sample = u_vertical_max_sample,
      .u_vertical_max_relative_error = relative_error(u_vertical_max_sample, reference.u_vertical_max),
      .u_vertical_min_sample = u_vertical_min_sample,
      .u_vertical_min_relative_error = relative_error(u_vertical_min_sample, reference.u_vertical_min),
      .v_horizontal_max_sample = v_horizontal_max_sample,
      .v_horizontal_max_relative_error =
          relative_error(v_horizontal_max_sample, reference.v_horizontal_max),
      .v_horizontal_min_sample = v_horizontal_min_sample,
      .v_horizontal_min_relative_error =
          relative_error(v_horizontal_min_sample, reference.v_horizontal_min),
      .divergence_l2 = result.final_step.divergence_l2,
  };
  validation.pass = validation.u_vertical_max_relative_error <= 0.02 &&
                    validation.u_vertical_min_relative_error <= 0.02 &&
                    validation.v_horizontal_max_relative_error <= 0.02 &&
                    validation.v_horizontal_min_relative_error <= 0.02 &&
                    validation.divergence_l2 <= 1.0e-10;
  return validation;
}

LidDrivenCavityResult run_lid_driven_cavity(const LidDrivenCavityConfig& config) {
  validate_config(config);

  const Grid grid{config.nx,
                  config.ny,
                  1,
                  1.0 / static_cast<double>(config.nx),
                  1.0 / static_cast<double>(config.ny),
                  1.0,
                  1};
  const BoundaryConditionSet boundary_conditions = make_lid_driven_cavity_boundary_conditions(config);
  const PressureBoundarySet pressure_boundary_conditions =
      derive_pressure_boundary_conditions(boundary_conditions);
  const double viscosity = config.lid_velocity / config.reynolds;
  const double dt = fixed_dt(config, grid);
  const double predictor_alpha = 0.5 * viscosity * dt;
  const ProjectionOptions projection_options{
      .dt = dt,
      .density = 1.0,
      .poisson_max_iterations = config.poisson_max_iterations,
      .poisson_tolerance = config.poisson_tolerance,
  };

  VelocityField velocity{grid};
  VelocityField advection_current{grid};
  VelocityField advection_previous{grid};
  VelocityField diffusion{grid};
  VelocityField pressure_gradient{grid};
  VelocityField predictor_rhs{grid};
  VelocityField predicted{grid};
  VelocityField corrected{grid};
  PressureField pressure_total{grid};
  PressureField pressure_correction{grid};
  ProjectionDiagnostics projection{};
  bool has_previous_advection = false;
  double time = 0.0;
  SimulationStepMetrics metrics{.dt = dt};

  velocity.fill(0.0);
  pressure_total.fill(0.0);
  pressure_correction.fill(0.0);
  apply_velocity_boundary_conditions(boundary_conditions, velocity);
  apply_pressure_boundary_conditions(pressure_boundary_conditions, pressure_total);

  for(int step = 0; step < config.max_steps; ++step) {
    VelocityField current = velocity;
    apply_velocity_boundary_conditions(boundary_conditions, current);

    compute_advection_term(current, config.advection, advection_current);
    compute_diffusion_term(current, viscosity, diffusion);
    operators::compute_gradient(pressure_total, pressure_gradient);

    predictor_rhs = current;
    axpy_velocity(predictor_rhs, pressure_gradient, -dt);
    axpy_velocity(predictor_rhs, diffusion, 0.5 * dt);
    if(has_previous_advection) {
      axpy_velocity(predictor_rhs, advection_current, -1.5 * dt);
      axpy_velocity(predictor_rhs, advection_previous, 0.5 * dt);
    } else {
      axpy_velocity(predictor_rhs, advection_current, -dt);
    }

    predicted = predictor_rhs;
    solve_predictor_adi(predictor_rhs, predictor_alpha, boundary_conditions, predicted);
    pressure_correction.fill(0.0);
    projection = project_velocity(
        predicted, boundary_conditions, projection_options, pressure_correction, corrected);
    axpy_active(pressure_total, pressure_correction, 1.0);
    apply_pressure_boundary_conditions(pressure_boundary_conditions, pressure_total);

    const CflDiagnostics cfl = compute_advective_cfl(corrected, dt);
    const double delta = max_velocity_change(velocity, corrected);
    if(!std::isfinite(delta) || !std::isfinite(projection.divergence_l2_after)) {
      throw std::runtime_error("lid-driven cavity solve produced a non-finite state");
    }

    velocity = corrected;
    advection_previous = advection_current;
    has_previous_advection = true;
    time += dt;

    metrics = SimulationStepMetrics{
        .step = step + 1,
        .time = time,
        .dt = dt,
        .max_cfl = cfl.max_cfl,
        .max_velocity_change = delta,
        .divergence_l2 = projection.divergence_l2_after,
        .pressure_iterations = projection.pressure_solve.iterations,
        .pressure_relative_residual = projection.pressure_solve.relative_residual,
    };

    if(metrics.step >= config.min_steps && metrics.max_velocity_change <= config.steady_tolerance) {
      break;
    }
  }

  apply_velocity_boundary_conditions(boundary_conditions, velocity);

  LidDrivenCavityResult result{
      .config = config,
      .final_step = metrics,
      .u_vertical_centerline = sample_vertical_centerline_u(velocity),
      .v_horizontal_centerline = sample_horizontal_centerline_v(velocity),
  };
  result.extrema = compute_centerline_extrema(result.u_vertical_centerline, result.v_horizontal_centerline);
  if(config.validate_reference && std::abs(config.reynolds - 100.0) < 1.0e-12) {
    result.validation = validate_lid_driven_cavity_re100(result);
  }
  return result;
}

}  // namespace solver
