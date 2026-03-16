#include "solver/taylor_green.hpp"

#include "operators/discrete_operators.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace solver {

namespace {

double pi() {
  return std::acos(-1.0);
}

double square(const double value) {
  return value * value;
}

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

void validate_config(const TaylorGreenConfig& config) {
  if(config.nx <= 0 || config.ny <= 0) {
    throw std::invalid_argument("Taylor-Green config expects positive nx and ny");
  }
  if(config.viscosity <= 0.0) {
    throw std::invalid_argument("Taylor-Green config expects viscosity > 0");
  }
  if(config.cfl_limit <= 0.0) {
    throw std::invalid_argument("Taylor-Green config expects cfl_limit > 0");
  }
  if(config.final_time <= 0.0) {
    throw std::invalid_argument("Taylor-Green config expects final_time > 0");
  }
  if(config.poisson_max_iterations <= 0 || config.poisson_tolerance <= 0.0) {
    throw std::invalid_argument("Taylor-Green config expects positive Poisson controls");
  }
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

double mixed_second_derivative_2d(const StructuredField& input, const int i, const int j, const int k) {
  const double dy = input.layout().grid().dy;
  const double dy_squared = dy * dy;
  const double dyy_im1 =
      (input(i - 1, j + 1, k) - 2.0 * input(i - 1, j, k) + input(i - 1, j - 1, k)) / dy_squared;
  const double dyy_i =
      (input(i, j + 1, k) - 2.0 * input(i, j, k) + input(i, j - 1, k)) / dy_squared;
  const double dyy_ip1 =
      (input(i + 1, j + 1, k) - 2.0 * input(i + 1, j, k) + input(i + 1, j - 1, k)) / dy_squared;

  const double dx = input.layout().grid().dx;
  const double inverse_dx_squared = 1.0 / (dx * dx);
  return (dyy_ip1 - 2.0 * dyy_i + dyy_im1) * inverse_dx_squared;
}

void compute_factorized_correction_2d(const VelocityField& current_velocity,
                                      VelocityField& factorized_correction) {
  factorized_correction.fill(0.0);

  const auto fill_component = [](const StructuredField& input, StructuredField& output) {
    const IndexRange3D active = input.layout().active_range();
    for(int k = active.k_begin; k < active.k_end; ++k) {
      for(int j = active.j_begin; j < active.j_end; ++j) {
        for(int i = active.i_begin; i < active.i_end; ++i) {
          output(i, j, k) = mixed_second_derivative_2d(input, i, j, k);
        }
      }
    }
  };

  fill_component(current_velocity.x, factorized_correction.x);
  fill_component(current_velocity.y, factorized_correction.y);
  fill_component(current_velocity.z, factorized_correction.z);
}

void assemble_predictor_rhs(const VelocityField& current_velocity,
                            const PressureField& pressure_total,
                            const VelocityField* previous_advection,
                            const AdvectionOptions& advection_options,
                            const double viscosity,
                            const double dt,
                            VelocityField& advection_current,
                            VelocityField& diffusion,
                            VelocityField& pressure_gradient,
                            VelocityField& factorized_correction,
                            VelocityField& predictor_rhs) {
  compute_advection_term(current_velocity, advection_options, advection_current);
  compute_diffusion_term(current_velocity, viscosity, diffusion);
  operators::compute_gradient(pressure_total, pressure_gradient);
  compute_factorized_correction_2d(current_velocity, factorized_correction);

  predictor_rhs = current_velocity;
  axpy_velocity(predictor_rhs, pressure_gradient, -dt);
  axpy_velocity(predictor_rhs, diffusion, 0.5 * dt);

  if(previous_advection != nullptr) {
    axpy_velocity(predictor_rhs, advection_current, -1.5 * dt);
    axpy_velocity(predictor_rhs, *previous_advection, 0.5 * dt);
  } else {
    axpy_velocity(predictor_rhs, advection_current, -dt);
  }

  const double alpha = 0.5 * viscosity * dt;
  axpy_velocity(predictor_rhs, factorized_correction, alpha * alpha);
}

double exact_decay(const TaylorGreenConfig& config, const double time) {
  return std::exp(-2.0 * config.viscosity * time);
}

double exact_kinetic_energy(const TaylorGreenConfig& config, const double time) {
  return 0.25 * std::exp(-4.0 * config.viscosity * time);
}

double kinetic_energy(const VelocityField& velocity) {
  const Grid& grid = velocity.x.layout().grid();
  const IndexRange3D cells = FieldLayout::cell_centered(grid).active_range();
  double energy_sum = 0.0;
  std::size_t count = 0;

  for(int k = cells.k_begin; k < cells.k_end; ++k) {
    for(int j = cells.j_begin; j < cells.j_end; ++j) {
      for(int i = cells.i_begin; i < cells.i_end; ++i) {
        const double u_center = 0.5 * (velocity.x(i, j, k) + velocity.x(i + 1, j, k));
        const double v_center = 0.5 * (velocity.y(i, j, k) + velocity.y(i, j + 1, k));
        const double w_center = 0.5 * (velocity.z(i, j, k) + velocity.z(i, j, k + 1));
        energy_sum += 0.5 * (square(u_center) + square(v_center) + square(w_center));
        ++count;
      }
    }
  }

  return energy_sum / static_cast<double>(count);
}

double velocity_relative_l2_error(const VelocityField& velocity,
                                  const TaylorGreenConfig& config,
                                  const double time) {
  const double decay = exact_decay(config, time);
  double numerator = 0.0;
  double denominator = 0.0;

  const auto accumulate = [&](const StructuredField& field, const auto& exact_value) {
    const IndexRange3D active = field.layout().active_range();
    for(int k = active.k_begin; k < active.k_end; ++k) {
      const double z = field.layout().coordinate_for_storage_index(Axis::z, k);
      for(int j = active.j_begin; j < active.j_end; ++j) {
        const double y = field.layout().coordinate_for_storage_index(Axis::y, j);
        for(int i = active.i_begin; i < active.i_end; ++i) {
          const double x = field.layout().coordinate_for_storage_index(Axis::x, i);
          const double reference = exact_value(x, y, z, decay);
          numerator += square(field(i, j, k) - reference);
          denominator += square(reference);
        }
      }
    }
  };

  accumulate(velocity.x, [](const double x, const double y, double, const double decay_value) {
    return -std::cos(x) * std::sin(y) * decay_value;
  });
  accumulate(velocity.y, [](const double x, const double y, double, const double decay_value) {
    return std::sin(x) * std::cos(y) * decay_value;
  });
  accumulate(velocity.z, [](double, double, double, double) {
    return 0.0;
  });

  return denominator > 0.0 ? std::sqrt(numerator / denominator) : std::sqrt(numerator);
}

void initialize_exact_state(const TaylorGreenConfig& config,
                            const BoundaryConditionSet& boundary_conditions,
                            VelocityField& velocity,
                            PressureField& pressure_total) {
  const double time = 0.0;
  const double decay = exact_decay(config, time);

  const IndexRange3D u_active = velocity.x.layout().active_range();
  for(int k = u_active.k_begin; k < u_active.k_end; ++k) {
    for(int j = u_active.j_begin; j < u_active.j_end; ++j) {
      const double y = velocity.x.layout().coordinate_for_storage_index(Axis::y, j);
      for(int i = u_active.i_begin; i < u_active.i_end; ++i) {
        const double x = velocity.x.layout().coordinate_for_storage_index(Axis::x, i);
        velocity.x(i, j, k) = -std::cos(x) * std::sin(y) * decay;
      }
    }
  }

  const IndexRange3D v_active = velocity.y.layout().active_range();
  for(int k = v_active.k_begin; k < v_active.k_end; ++k) {
    for(int j = v_active.j_begin; j < v_active.j_end; ++j) {
      const double y = velocity.y.layout().coordinate_for_storage_index(Axis::y, j);
      for(int i = v_active.i_begin; i < v_active.i_end; ++i) {
        const double x = velocity.y.layout().coordinate_for_storage_index(Axis::x, i);
        velocity.y(i, j, k) = std::sin(x) * std::cos(y) * decay;
      }
    }
  }

  velocity.z.fill(0.0);

  const IndexRange3D p_active = pressure_total.layout().active_range();
  const double pressure_decay = std::exp(-4.0 * config.viscosity * time);
  for(int k = p_active.k_begin; k < p_active.k_end; ++k) {
    for(int j = p_active.j_begin; j < p_active.j_end; ++j) {
      const double y = pressure_total.layout().coordinate_for_storage_index(Axis::y, j);
      for(int i = p_active.i_begin; i < p_active.i_end; ++i) {
        const double x = pressure_total.layout().coordinate_for_storage_index(Axis::x, i);
        pressure_total(i, j, k) = -0.25 * (std::cos(2.0 * x) + std::cos(2.0 * y)) * pressure_decay;
      }
    }
  }

  apply_velocity_boundary_conditions(boundary_conditions, velocity);
  VelocityField diffusion{velocity.x.layout().grid()};
  compute_diffusion_term(velocity, config.viscosity, diffusion);
  apply_total_pressure_boundary_conditions(boundary_conditions, diffusion, pressure_total);
}

double stable_dt_upper_bound(const Grid& grid, const TaylorGreenConfig& config) {
  return config.cfl_limit / (1.0 / grid.dx + 1.0 / grid.dy);
}

TaylorGreenValidation build_validation(const TaylorGreenConfig& config,
                                       const TaylorGreenStepMetrics& metrics,
                                       const VelocityField& velocity) {
  const double exact_energy = exact_kinetic_energy(config, metrics.time);
  const double energy_error =
      std::abs(kinetic_energy(velocity) - exact_energy) / exact_energy;
  return TaylorGreenValidation{
      .reference_dataset = "analytic_taylor_green_decay",
      .normalized_energy_error = energy_error,
      .velocity_relative_l2_error = velocity_relative_l2_error(velocity, config, metrics.time),
      .exact_kinetic_energy = exact_energy,
      .pass = energy_error <= 1.0e-2 && metrics.divergence_l2 <= 1.0e-10,
  };
}

}  // namespace

TaylorGreenConfig default_taylor_green_config() {
  return TaylorGreenConfig{};
}

TaylorGreenConfig load_taylor_green_config(const std::string& path) {
  std::ifstream input(path);
  if(!input.is_open()) {
    throw std::runtime_error("unable to open Taylor-Green config: " + path);
  }

  TaylorGreenConfig config = default_taylor_green_config();
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
    } else if(key == "viscosity") {
      config.viscosity = std::stod(value);
    } else if(key == "cfl_limit") {
      config.cfl_limit = std::stod(value);
    } else if(key == "final_time") {
      config.final_time = std::stod(value);
    } else if(key == "poisson_max_iterations") {
      config.poisson_max_iterations = std::stoi(value);
    } else if(key == "poisson_tolerance") {
      config.poisson_tolerance = std::stod(value);
    } else if(key == "validate_energy") {
      config.validate_energy = parse_bool(value);
    } else {
      throw std::runtime_error("unsupported Taylor-Green config key: " + key);
    }
  }

  validate_config(config);
  return config;
}

std::string describe(const TaylorGreenConfig& config) {
  std::ostringstream builder;
  builder << "nx=" << config.nx << ", ny=" << config.ny
          << ", viscosity=" << config.viscosity
          << ", cfl_limit=" << config.cfl_limit
          << ", final_time=" << config.final_time
          << ", poisson_max_iterations=" << config.poisson_max_iterations
          << ", poisson_tolerance=" << config.poisson_tolerance
          << ", validate_energy=" << (config.validate_energy ? "true" : "false")
          << ", advection=" << describe(config.advection);
  return builder.str();
}

BoundaryConditionSet make_taylor_green_boundary_conditions() {
  BoundaryConditionSet boundary_conditions = BoundaryConditionSet::all(PhysicalBoundaryType::periodic);
  boundary_conditions[BoundaryFace::z_min].type = PhysicalBoundaryType::symmetry;
  boundary_conditions[BoundaryFace::z_max].type = PhysicalBoundaryType::symmetry;
  return boundary_conditions;
}

TaylorGreenState::TaylorGreenState(const Grid& grid_in)
    : grid(grid_in),
      velocity(grid_in),
      advection_previous(grid_in),
      pressure_total(grid_in) {}

double taylor_green_dt(const TaylorGreenConfig& config) {
  validate_config(config);

  const double domain_length = 2.0 * pi();
  const Grid grid{config.nx,
                  config.ny,
                  1,
                  domain_length / static_cast<double>(config.nx),
                  domain_length / static_cast<double>(config.ny),
                  1.0,
                  1};
  const double dt_upper = stable_dt_upper_bound(grid, config);
  const int step_count = std::max(1, static_cast<int>(std::ceil(config.final_time / dt_upper)));
  return config.final_time / static_cast<double>(step_count);
}

TaylorGreenState initialize_taylor_green_state(const TaylorGreenConfig& config) {
  validate_config(config);
  const double domain_length = 2.0 * pi();
  const Grid grid{config.nx,
                  config.ny,
                  1,
                  domain_length / static_cast<double>(config.nx),
                  domain_length / static_cast<double>(config.ny),
                  1.0,
                  1};
  TaylorGreenState state{grid};
  state.metrics.dt = taylor_green_dt(config);
  initialize_exact_state(config, make_taylor_green_boundary_conditions(), state.velocity, state.pressure_total);
  state.advection_previous.fill(0.0);
  return state;
}

void run_taylor_green_steps(const TaylorGreenConfig& config,
                            const int step_count,
                            TaylorGreenState& state) {
  validate_config(config);
  if(step_count < 0) {
    throw std::invalid_argument("run_taylor_green_steps expects a non-negative step count");
  }

  const double domain_length = 2.0 * pi();
  const Grid expected_grid{config.nx,
                           config.ny,
                           1,
                           domain_length / static_cast<double>(config.nx),
                           domain_length / static_cast<double>(config.ny),
                           1.0,
                           1};
  if(state.grid.nx != expected_grid.nx || state.grid.ny != expected_grid.ny ||
     std::abs(state.grid.dx - expected_grid.dx) > 0.0 || std::abs(state.grid.dy - expected_grid.dy) > 0.0) {
    throw std::invalid_argument("Taylor-Green state grid does not match the config");
  }

  const BoundaryConditionSet boundary_conditions = make_taylor_green_boundary_conditions();
  const double dt = taylor_green_dt(config);
  const double predictor_alpha = 0.5 * config.viscosity * dt;
  const ProjectionOptions projection_options{
      .dt = dt,
      .density = 1.0,
      .poisson_max_iterations = config.poisson_max_iterations,
      .poisson_tolerance = config.poisson_tolerance,
  };

  VelocityField advection_current{state.grid};
  VelocityField diffusion{state.grid};
  VelocityField pressure_gradient{state.grid};
  VelocityField factorized_correction{state.grid};
  VelocityField predictor_rhs{state.grid};
  VelocityField predicted{state.grid};
  VelocityField corrected{state.grid};
  PressureField pressure_correction{state.grid};
  ScalarField pressure_rhs{state.grid};
  ProjectionDiagnostics projection{};

  for(int step = 0; step < step_count; ++step) {
    VelocityField current = state.velocity;
    apply_velocity_boundary_conditions(boundary_conditions, current);

    assemble_predictor_rhs(current,
                           state.pressure_total,
                           state.has_previous_advection ? &state.advection_previous : nullptr,
                           config.advection,
                           config.viscosity,
                           dt,
                           advection_current,
                           diffusion,
                           pressure_gradient,
                           factorized_correction,
                           predictor_rhs);

    predicted = predictor_rhs;
    solve_predictor_adi(predictor_rhs, predictor_alpha, boundary_conditions, predicted);

    pressure_correction.fill(0.0);
    projection = project_velocity(
        predicted,
        boundary_conditions,
        projection_options,
        pressure_correction,
        corrected,
        &pressure_rhs);

    axpy_active(state.pressure_total, pressure_correction, 1.0);
    axpy_active(state.pressure_total, pressure_rhs, -0.5 * config.viscosity * dt);
    compute_diffusion_term(corrected, config.viscosity, diffusion);
    apply_total_pressure_boundary_conditions(boundary_conditions, diffusion, state.pressure_total);

    const CflDiagnostics cfl = compute_advective_cfl(corrected, dt);
    const double delta = max_velocity_change(state.velocity, corrected);
    if(!std::isfinite(delta) || !std::isfinite(projection.divergence_l2_after)) {
      throw std::runtime_error("Taylor-Green solve produced a non-finite state");
    }

    state.velocity = corrected;
    state.advection_previous = advection_current;
    state.has_previous_advection = true;

    state.metrics = TaylorGreenStepMetrics{
        .step = state.metrics.step + 1,
        .time = state.metrics.time + dt,
        .dt = dt,
        .max_cfl = cfl.max_cfl,
        .max_velocity_change = delta,
        .divergence_l2 = projection.divergence_l2_after,
        .max_divergence_l2 = std::max(state.metrics.max_divergence_l2, projection.divergence_l2_after),
        .pressure_iterations = projection.pressure_solve.iterations,
        .pressure_relative_residual = projection.pressure_solve.relative_residual,
    };
  }
}

TaylorGreenResult finalize_taylor_green_result(const TaylorGreenConfig& config,
                                               const TaylorGreenState& state) {
  validate_config(config);
  TaylorGreenResult result{
      .config = config,
      .final_step = state.metrics,
      .initial_kinetic_energy = exact_kinetic_energy(config, 0.0),
      .final_kinetic_energy = kinetic_energy(state.velocity),
  };
  if(config.validate_energy) {
    result.validation = build_validation(config, state.metrics, state.velocity);
  }
  return result;
}

TaylorGreenResult run_taylor_green(const TaylorGreenConfig& config) {
  TaylorGreenState state = initialize_taylor_green_state(config);
  run_taylor_green_steps(config, static_cast<int>(std::ceil(config.final_time / taylor_green_dt(config))), state);
  return finalize_taylor_green_result(config, state);
}

}  // namespace solver
