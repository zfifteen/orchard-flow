#include "core/fields.hpp"
#include "core/grid.hpp"
#include "core/runtime.hpp"
#include "linsolve/poisson_solver.hpp"
#include "operators/discrete_operators.hpp"
#include "solver/lid_driven_cavity.hpp"
#include "solver/momentum_terms.hpp"
#include "solver/projection.hpp"

#include <cmath>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>
#include <type_traits>

namespace {

void require(bool condition, const std::string& message) {
  if(!condition) {
    throw std::runtime_error(message);
  }
}

double pi() {
  return std::acos(-1.0);
}

std::string source_path(const std::string& relative_path) {
  return std::string(SOLVER_SOURCE_DIR) + "/" + relative_path;
}

double square(const double value) {
  return value * value;
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

template <typename Field, typename ValueFn>
double active_l2_error(const Field& field, ValueFn&& exact_value) {
  const solver::IndexRange3D active = field.layout().active_range();

  double sum = 0.0;
  std::size_t count = 0;
  for(int k = active.k_begin; k < active.k_end; ++k) {
    const double z = field.layout().coordinate_for_storage_index(solver::Axis::z, k);
    for(int j = active.j_begin; j < active.j_end; ++j) {
      const double y = field.layout().coordinate_for_storage_index(solver::Axis::y, j);
      for(int i = active.i_begin; i < active.i_end; ++i) {
        const double x = field.layout().coordinate_for_storage_index(solver::Axis::x, i);
        sum += square(field(i, j, k) - exact_value(x, y, z));
        ++count;
      }
    }
  }

  return std::sqrt(sum / static_cast<double>(count));
}

double observed_order(const double coarse_error, const double fine_error) {
  return std::log(coarse_error / fine_error) / std::log(2.0);
}

double wrap_periodic(const double coordinate, const double period) {
  const double wrapped = std::fmod(coordinate, period);
  return wrapped < 0.0 ? wrapped + period : wrapped;
}

double kinetic_energy(const solver::VelocityField& velocity) {
  const solver::Grid& grid = velocity.x.layout().grid();
  const solver::IndexRange3D cells = solver::FieldLayout::cell_centered(grid).active_range();
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

double active_min(const solver::StructuredField& field) {
  const solver::IndexRange3D active = field.layout().active_range();
  double value = field(active.i_begin, active.j_begin, active.k_begin);

  for(int k = active.k_begin; k < active.k_end; ++k) {
    for(int j = active.j_begin; j < active.j_end; ++j) {
      for(int i = active.i_begin; i < active.i_end; ++i) {
        value = std::min(value, field(i, j, k));
      }
    }
  }

  return value;
}

double active_max(const solver::StructuredField& field) {
  const solver::IndexRange3D active = field.layout().active_range();
  double value = field(active.i_begin, active.j_begin, active.k_begin);

  for(int k = active.k_begin; k < active.k_end; ++k) {
    for(int j = active.j_begin; j < active.j_end; ++j) {
      for(int i = active.i_begin; i < active.i_end; ++i) {
        value = std::max(value, field(i, j, k));
      }
    }
  }

  return value;
}

double max_abs_active(const solver::StructuredField& field) {
  const solver::IndexRange3D active = field.layout().active_range();
  double value = 0.0;

  for(int k = active.k_begin; k < active.k_end; ++k) {
    for(int j = active.j_begin; j < active.j_end; ++j) {
      for(int i = active.i_begin; i < active.i_end; ++i) {
        value = std::max(value, std::abs(field(i, j, k)));
      }
    }
  }

  return value;
}

double active_mean_value(const solver::StructuredField& field) {
  const solver::IndexRange3D active = field.layout().active_range();
  const std::size_t count = active.extent().cell_count();
  double sum = 0.0;

  for(int k = active.k_begin; k < active.k_end; ++k) {
    for(int j = active.j_begin; j < active.j_end; ++j) {
      for(int i = active.i_begin; i < active.i_end; ++i) {
        sum += field(i, j, k);
      }
    }
  }

  return sum / static_cast<double>(count);
}

template <typename Field>
double active_l2_difference(const Field& left, const Field& right) {
  const solver::IndexRange3D active = left.layout().active_range();
  double sum = 0.0;
  std::size_t count = 0;

  for(int k = active.k_begin; k < active.k_end; ++k) {
    for(int j = active.j_begin; j < active.j_end; ++j) {
      for(int i = active.i_begin; i < active.i_end; ++i) {
        sum += square(left(i, j, k) - right(i, j, k));
        ++count;
      }
    }
  }

  return std::sqrt(sum / static_cast<double>(count));
}

template <typename Field>
double active_l2_norm_value(const Field& field) {
  const solver::IndexRange3D active = field.layout().active_range();
  double sum = 0.0;
  std::size_t count = 0;

  for(int k = active.k_begin; k < active.k_end; ++k) {
    for(int j = active.j_begin; j < active.j_end; ++j) {
      for(int i = active.i_begin; i < active.i_end; ++i) {
        sum += square(field(i, j, k));
        ++count;
      }
    }
  }

  return std::sqrt(sum / static_cast<double>(count));
}

template <typename Field>
double relative_active_l2_difference(const Field& left, const Field& right) {
  const double denominator = active_l2_norm_value(right);
  require(denominator > 0.0, "relative_active_l2_difference requires non-zero reference norm");
  return active_l2_difference(left, right) / denominator;
}

template <typename Field>
void subtract_active_mean_inplace(Field& field) {
  const double mean = active_mean_value(field);
  const solver::IndexRange3D active = field.layout().active_range();

  for(int k = active.k_begin; k < active.k_end; ++k) {
    for(int j = active.j_begin; j < active.j_end; ++j) {
      for(int i = active.i_begin; i < active.i_end; ++i) {
        field(i, j, k) -= mean;
      }
    }
  }
}

template <typename FromField, typename ToField>
void copy_active_values(const FromField& source, ToField& destination) {
  const solver::IndexRange3D active = source.layout().active_range();

  for(int k = active.k_begin; k < active.k_end; ++k) {
    for(int j = active.j_begin; j < active.j_end; ++j) {
      for(int i = active.i_begin; i < active.i_end; ++i) {
        destination(i, j, k) = source(i, j, k);
      }
    }
  }
}

template <typename Field>
void axpy_active(Field& destination, const Field& source, const double scale) {
  const solver::IndexRange3D active = destination.layout().active_range();

  for(int k = active.k_begin; k < active.k_end; ++k) {
    for(int j = active.j_begin; j < active.j_end; ++j) {
      for(int i = active.i_begin; i < active.i_end; ++i) {
        destination(i, j, k) += scale * source(i, j, k);
      }
    }
  }
}

void axpy_velocity(solver::VelocityField& destination,
                   const solver::VelocityField& source,
                   const double scale) {
  axpy_active(destination.x, source.x, scale);
  axpy_active(destination.y, source.y, scale);
  axpy_active(destination.z, source.z, scale);
}

template <typename Field>
void scale_active(Field& field, const double scale) {
  const solver::IndexRange3D active = field.layout().active_range();

  for(int k = active.k_begin; k < active.k_end; ++k) {
    for(int j = active.j_begin; j < active.j_end; ++j) {
      for(int i = active.i_begin; i < active.i_end; ++i) {
        field(i, j, k) *= scale;
      }
    }
  }
}

void scale_velocity_active(solver::VelocityField& field, const double scale) {
  scale_active(field.x, scale);
  scale_active(field.y, scale);
  scale_active(field.z, scale);
}

double velocity_max_abs(const solver::VelocityField& velocity) {
  return std::max({max_abs_active(velocity.x), max_abs_active(velocity.y), max_abs_active(velocity.z)});
}

void test_build_profile_is_locked() {
  const solver::BuildInfo build_info = solver::get_build_info();
  require(build_info.build_profile == "deterministic" || build_info.build_profile == "benchmark",
          "unexpected build profile");
}

void test_runtime_platform_is_supported() {
  const solver::BuildInfo build_info = solver::get_build_info();
  require(build_info.supported_runtime_platform, "expected Apple Silicon runtime");
}

void test_banner_contains_profile() {
  const solver::BuildInfo build_info = solver::get_build_info();
  const std::string banner = solver::format_build_banner(build_info);
  require(banner.find("profile: " + build_info.build_profile) != std::string::npos,
          "banner missing build profile");
}

void test_grid_coordinates() {
  const solver::Grid grid{4, 3, 2, 0.5, 0.25, 0.125, 1};

  require(!grid.is_2d(), "expected 3D grid");
  require(grid.cells(solver::Axis::x) == 4, "wrong x cell count");
  require(grid.cell_center(solver::Axis::x, 0) == 0.25, "wrong x cell center");
  require(grid.cell_center(solver::Axis::y, 1) == 0.375, "wrong y cell center");
  require(grid.face_coordinate(solver::Axis::z, 2) == 0.25, "wrong z face coordinate");
}

void test_pressure_layout_indexing() {
  const solver::Grid grid{4, 3, 2, 0.5, 0.25, 0.125, 1};
  const solver::PressureField pressure{grid};
  const solver::FieldLayout& layout = pressure.layout();
  const solver::Extent3D active = layout.active_extent();
  const solver::Extent3D storage = layout.storage_extent();
  const solver::IndexRange3D active_range = layout.active_range();

  require(active.nx == 4 && active.ny == 3 && active.nz == 2, "wrong pressure active extent");
  require(storage.nx == 6 && storage.ny == 5 && storage.nz == 4, "wrong pressure storage extent");
  require(active_range.i_begin == 1 && active_range.i_end == 5, "wrong active i range");
  require(active_range.j_begin == 1 && active_range.j_end == 4, "wrong active j range");
  require(active_range.k_begin == 1 && active_range.k_end == 3, "wrong active k range");
  require(layout.is_unit_stride_i(), "expected i to be unit stride");

  const std::size_t expected =
      static_cast<std::size_t>(2) + static_cast<std::size_t>(storage.nx) *
                                       (static_cast<std::size_t>(3) +
                                        static_cast<std::size_t>(storage.ny) *
                                            static_cast<std::size_t>(1));
  require(layout.index(2, 3, 1) == expected, "flat index mapping does not match contract");
}

void test_ghost_cell_access_and_boundary_ranges() {
  const solver::Grid grid{4, 3, 2, 0.5, 0.25, 0.125, 1};
  solver::PressureField pressure{grid};
  pressure.fill(0.0);
  pressure.fill_ghost_layer(solver::BoundaryFace::x_min, 0, 42.0);

  const solver::IndexRange3D x_min_ghost = pressure.layout().ghost_range(solver::BoundaryFace::x_min);
  const solver::IndexRange3D x_max_boundary =
      pressure.layout().boundary_active_range(solver::BoundaryFace::x_max);

  for(int k = x_min_ghost.k_begin; k < x_min_ghost.k_end; ++k) {
    for(int j = x_min_ghost.j_begin; j < x_min_ghost.j_end; ++j) {
      for(int i = x_min_ghost.i_begin; i < x_min_ghost.i_end; ++i) {
        require(pressure(i, j, k) == 42.0, "ghost fill helper missed a cell");
      }
    }
  }

  require(x_max_boundary.i_begin == pressure.layout().active_range().i_end - 1,
          "wrong x_max boundary slab");
  require(x_max_boundary.i_end == pressure.layout().active_range().i_end,
          "wrong x_max boundary slab end");

  const solver::Index3D first_active = pressure.layout().storage_index_from_active(0, 0, 0);
  require(pressure(first_active.i, first_active.j, first_active.k) == 0.0,
          "ghost fill touched active storage");
}

void test_memory_layout_and_alignment() {
  const solver::Grid grid{4, 3, 2, 0.5, 0.25, 0.125, 1};
  solver::PressureField pressure{grid};
  const solver::IndexRange3D active_range = pressure.layout().active_range();

  require(pressure.is_aligned(), "pressure storage is not aligned");
  require(&pressure(active_range.i_begin + 1, active_range.j_begin, active_range.k_begin) ==
              &pressure(active_range.i_begin, active_range.j_begin, active_range.k_begin) + 1,
          "i-direction is not contiguous");
}

void test_cell_and_face_placement() {
  const solver::Grid grid{4, 3, 2, 0.5, 0.25, 0.125, 1};
  const solver::PressureField pressure{grid};
  const solver::VelocityField velocity{grid};

  require(pressure.layout().active_extent().nx == 4, "wrong pressure x extent");
  require(velocity.x.layout().active_extent().nx == 5, "wrong u extent");
  require(velocity.y.layout().active_extent().ny == 4, "wrong v extent");
  require(velocity.z.layout().active_extent().nz == 3, "wrong w extent");

  require(pressure.layout().coordinate_at_active_index(solver::Axis::x, 0) == 0.25,
          "wrong pressure x coordinate");
  require(velocity.x.layout().coordinate_at_active_index(solver::Axis::x, 0) == 0.0,
          "wrong u-face x coordinate");
  require(velocity.x.layout().coordinate_at_active_index(solver::Axis::y, 0) == 0.125,
          "wrong u-face y coordinate");
  require(velocity.y.layout().coordinate_at_active_index(solver::Axis::y, 0) == 0.0,
          "wrong v-face y coordinate");
  require(velocity.z.layout().coordinate_at_active_index(solver::Axis::z, 0) == 0.0,
          "wrong w-face z coordinate");
}

void test_double_precision_storage_contract() {
  static_assert(std::is_same_v<solver::StructuredField::value_type, double>,
                "structured fields must store doubles");
  static_assert(std::is_same_v<solver::PressureField::value_type, double>,
                "pressure field must store doubles");
  static_assert(std::is_same_v<solver::ScalarField::value_type, double>,
                "scalar field must store doubles");

  const solver::Grid grid{4, 3, 2, 0.5, 0.25, 0.125, 1};
  const solver::PressureField pressure{grid};

  require(sizeof(*pressure.data()) == sizeof(double), "unexpected pressure storage precision");
}

void test_advection_options_and_cfl_diagnostic() {
  const solver::AdvectionOptions options{};
  require(solver::to_string(options.scheme) == "tvd", "wrong default advection scheme label");
  require(solver::to_string(options.limiter) == "van_leer", "wrong default limiter label");
  require(solver::describe(options).find("scheme=tvd") != std::string::npos,
          "missing scheme description");
  require(solver::describe(options).find("limiter=van_leer") != std::string::npos,
          "missing limiter description");

  const solver::Grid grid{8, 4, 1, 0.25, 0.5, 1.0, 1};
  solver::VelocityField velocity{grid};
  fill_storage(velocity.x, [](double, double, double) { return 2.0; });
  fill_storage(velocity.y, [](double, double, double) { return -1.0; });
  fill_storage(velocity.z, [](double, double, double) { return 0.0; });

  const double dt = 0.05;
  const solver::CflDiagnostics diagnostics = solver::compute_advective_cfl(velocity, dt);
  require(std::abs(diagnostics.max_u - 2.0) < 1.0e-12, "wrong max u in CFL diagnostic");
  require(std::abs(diagnostics.max_v - 1.0) < 1.0e-12, "wrong max v in CFL diagnostic");
  require(std::abs(diagnostics.max_w) < 1.0e-12, "wrong max w in CFL diagnostic");
  require(std::abs(diagnostics.max_cfl - dt * (2.0 / 0.25 + 1.0 / 0.5)) < 1.0e-12,
          "wrong CFL value");
}

void test_projection_boundary_mapping() {
  solver::BoundaryConditionSet boundary_conditions{};
  boundary_conditions[solver::BoundaryFace::x_min].type = solver::PhysicalBoundaryType::no_slip_wall;
  boundary_conditions[solver::BoundaryFace::x_max].type =
      solver::PhysicalBoundaryType::prescribed_velocity;
  boundary_conditions[solver::BoundaryFace::y_min].type = solver::PhysicalBoundaryType::symmetry;
  boundary_conditions[solver::BoundaryFace::y_max].type = solver::PhysicalBoundaryType::fixed_pressure;
  boundary_conditions[solver::BoundaryFace::y_max].pressure = 1.25;
  boundary_conditions[solver::BoundaryFace::z_min].type = solver::PhysicalBoundaryType::periodic;
  boundary_conditions[solver::BoundaryFace::z_max].type = solver::PhysicalBoundaryType::periodic;

  const solver::PressureBoundarySet mapped =
      solver::derive_pressure_boundary_conditions(boundary_conditions);

  require(mapped[solver::BoundaryFace::x_min].type == solver::PressureBoundaryType::neumann,
          "no-slip wall should map to pressure Neumann");
  require(mapped[solver::BoundaryFace::x_max].type == solver::PressureBoundaryType::neumann,
          "prescribed velocity should map to pressure Neumann");
  require(mapped[solver::BoundaryFace::y_min].type == solver::PressureBoundaryType::neumann,
          "symmetry should map to pressure Neumann");
  require(mapped[solver::BoundaryFace::y_max].type == solver::PressureBoundaryType::dirichlet,
          "fixed pressure should map to pressure Dirichlet");
  require(std::abs(mapped[solver::BoundaryFace::y_max].value - 1.25) < 1.0e-12,
          "fixed pressure value did not carry through the mapping");
  require(mapped[solver::BoundaryFace::z_min].type == solver::PressureBoundaryType::periodic,
          "periodic boundary should map to periodic pressure");
  require(mapped[solver::BoundaryFace::z_max].type == solver::PressureBoundaryType::periodic,
          "periodic boundary should map to periodic pressure");
}

void test_predictor_adi_preserves_quiescent_state() {
  const solver::Grid grid{12, 10, 1, 1.0 / 12.0, 1.0 / 10.0, 1.0, 1};
  const solver::BoundaryConditionSet boundary_conditions = solver::BoundaryConditionSet::cavity();

  solver::VelocityField rhs{grid};
  solver::VelocityField predicted{grid};
  rhs.fill(0.0);
  predicted.fill(1.0);

  const solver::HelmholtzDiagnostics diagnostics =
      solver::solve_predictor_adi(rhs, 0.05, boundary_conditions, predicted);

  require(diagnostics.line_solves > 0, "ADI predictor should perform deterministic line solves");
  require(velocity_max_abs(predicted) <= 1.0e-12, "quiescent predictor state should remain zero");
}

void test_poisson_mgpcg_discrete_dirichlet_recovery() {
  const int resolution = 32;
  const solver::Grid grid{resolution,
                          resolution,
                          1,
                          1.0 / static_cast<double>(resolution),
                          1.0 / static_cast<double>(resolution),
                          1.0,
                          1};
  solver::PressureBoundarySet boundary_conditions{};
  boundary_conditions[solver::BoundaryFace::x_min].type = solver::PressureBoundaryType::dirichlet;
  boundary_conditions[solver::BoundaryFace::x_max].type = solver::PressureBoundaryType::dirichlet;
  boundary_conditions[solver::BoundaryFace::y_min].type = solver::PressureBoundaryType::dirichlet;
  boundary_conditions[solver::BoundaryFace::y_max].type = solver::PressureBoundaryType::dirichlet;
  boundary_conditions[solver::BoundaryFace::z_min].type = solver::PressureBoundaryType::neumann;
  boundary_conditions[solver::BoundaryFace::z_max].type = solver::PressureBoundaryType::neumann;

  solver::PressureField exact_pressure{grid};
  solver::ScalarField rhs{grid};
  solver::PressureField solution{grid};
  fill_storage(exact_pressure, [](const double x, const double y, double) {
    return std::sin(pi() * x) * std::sin(pi() * y);
  });

  solver::linsolve::build_poisson_rhs_from_pressure(exact_pressure, boundary_conditions, rhs);
  solution.fill(0.0);

  const solver::ProjectionOptions options{
      .dt = 1.0,
      .density = 1.0,
      .poisson_max_iterations = 200,
      .poisson_tolerance = 1.0e-10,
  };
  const solver::PoissonSolveDiagnostics diagnostics =
      solver::linsolve::solve_pressure_poisson(rhs, boundary_conditions, options, solution);

  const double relative_error = relative_active_l2_difference(solution, exact_pressure);
  require(diagnostics.converged, "MGPCG Dirichlet solve should converge");
  require(diagnostics.relative_residual <= 1.0e-10, "MGPCG relative residual is too large");
  require(relative_error <= 1.0e-8, "MGPCG Dirichlet recovery error is too large");
  require(diagnostics.multigrid_levels >= 2, "MGPCG should build a multigrid hierarchy");
  require(diagnostics.coarse_unknowns > 0, "MGPCG should report the coarse-grid size");
  require(diagnostics.solver == "mgpcg", "unexpected pressure solver label");
  require(diagnostics.preconditioner == "geometric_multigrid", "unexpected preconditioner label");
  require(diagnostics.cycle == "v_cycle", "unexpected multigrid cycle label");
  require(diagnostics.smoother == "damped_jacobi", "unexpected smoother label");
  require(diagnostics.pre_smoothing_steps == 2, "unexpected pre-smoothing policy");
  require(diagnostics.post_smoothing_steps == 2, "unexpected post-smoothing policy");
}

void test_poisson_mgpcg_pure_neumann_zero_mean_recovery() {
  const int resolution = 24;
  const solver::Grid grid{resolution,
                          resolution,
                          1,
                          1.0 / static_cast<double>(resolution),
                          1.0 / static_cast<double>(resolution),
                          1.0,
                          1};
  solver::PressureBoundarySet boundary_conditions{};
  for(solver::PressureBoundaryCondition& boundary : boundary_conditions.faces) {
    boundary.type = solver::PressureBoundaryType::neumann;
  }

  solver::PressureField exact_pressure{grid};
  solver::ScalarField rhs{grid};
  solver::PressureField solution{grid};
  fill_storage(exact_pressure, [](const double x, const double y, double) {
    return std::cos(2.0 * pi() * x) * std::cos(2.0 * pi() * y);
  });
  subtract_active_mean_inplace(exact_pressure);

  solver::linsolve::build_poisson_rhs_from_pressure(exact_pressure, boundary_conditions, rhs);
  solution.fill(0.0);

  const solver::ProjectionOptions options{
      .dt = 1.0,
      .density = 1.0,
      .poisson_max_iterations = 300,
      .poisson_tolerance = 1.0e-10,
  };
  const solver::PoissonSolveDiagnostics diagnostics =
      solver::linsolve::solve_pressure_poisson(rhs, boundary_conditions, options, solution);

  const double relative_error = relative_active_l2_difference(solution, exact_pressure);
  require(diagnostics.converged, "MGPCG pure-Neumann solve should converge");
  require(diagnostics.zero_mean_enforced, "pure-Neumann solve should enforce zero mean");
  require(std::abs(active_mean_value(solution)) <= 1.0e-10,
          "pure-Neumann solve should preserve zero-mean pressure");
  require(diagnostics.relative_residual <= 1.0e-10, "pure-Neumann relative residual is too large");
  require(relative_error <= 1.0e-8, "pure-Neumann recovery error is too large");
}

struct ManufacturedErrors {
  double gradient;
  double divergence;
  double laplacian;
};

ManufacturedErrors run_manufactured_solution_case(const int resolution) {
  const double domain_length = 2.0 * pi();
  const double spacing = domain_length / static_cast<double>(resolution);
  const solver::Grid grid{resolution, resolution, 1, spacing, spacing, 1.0, 1};

  solver::PressureField pressure{grid};
  solver::VelocityField exact_velocity{grid};
  solver::VelocityField gradient{grid};
  solver::ScalarField divergence{grid};
  solver::PressureField laplacian{grid};

  fill_storage(pressure, [](const double x, const double y, double) {
    return std::sin(x) * std::cos(y);
  });

  fill_storage(exact_velocity.x, [](const double x, const double y, double) {
    return std::sin(x) * std::cos(y);
  });
  fill_storage(exact_velocity.y, [](const double x, const double y, double) {
    return std::cos(x) * std::sin(y);
  });
  fill_storage(exact_velocity.z, [](double, double, double) {
    return 0.0;
  });

  solver::operators::compute_gradient(pressure, gradient);
  solver::operators::compute_divergence(exact_velocity, divergence);
  solver::operators::compute_laplacian(pressure, laplacian);

  const double gradient_x_error = active_l2_error(gradient.x, [](const double x, const double y, double) {
    return std::cos(x) * std::cos(y);
  });
  const double gradient_y_error = active_l2_error(gradient.y, [](const double x, const double y, double) {
    return -std::sin(x) * std::sin(y);
  });
  const double divergence_error = active_l2_error(divergence, [](const double x, const double y, double) {
    return 2.0 * std::cos(x) * std::cos(y);
  });
  const double laplacian_error = active_l2_error(laplacian, [](const double x, const double y, double) {
    return -2.0 * std::sin(x) * std::cos(y);
  });

  return ManufacturedErrors{
      .gradient = std::sqrt(0.5 * (square(gradient_x_error) + square(gradient_y_error))),
      .divergence = divergence_error,
      .laplacian = laplacian_error,
  };
}

void test_manufactured_solution_convergence() {
  const ManufacturedErrors coarse = run_manufactured_solution_case(16);
  const ManufacturedErrors fine = run_manufactured_solution_case(32);

  require(coarse.gradient > fine.gradient, "gradient error did not decrease on refinement");
  require(coarse.divergence > fine.divergence, "divergence error did not decrease on refinement");
  require(coarse.laplacian > fine.laplacian, "laplacian error did not decrease on refinement");

  require(observed_order(coarse.gradient, fine.gradient) >= 1.8,
          "gradient convergence order below second-order target");
  require(observed_order(coarse.divergence, fine.divergence) >= 1.8,
          "divergence convergence order below second-order target");
  require(observed_order(coarse.laplacian, fine.laplacian) >= 1.8,
          "laplacian convergence order below second-order target");
}

void test_diffusion_term_matches_scaled_laplacian() {
  const double domain_length = 2.0 * pi();
  const int resolution = 32;
  const double spacing = domain_length / static_cast<double>(resolution);
  const solver::Grid grid{resolution, resolution, 1, spacing, spacing, 1.0, 1};
  const double viscosity = 0.125;

  solver::VelocityField velocity{grid};
  solver::VelocityField diffusion{grid};

  fill_storage(velocity.x, [](const double x, const double y, double) {
    return -std::cos(x) * std::sin(y);
  });
  fill_storage(velocity.y, [](const double x, const double y, double) {
    return std::sin(x) * std::cos(y);
  });
  fill_storage(velocity.z, [](double, double, double) {
    return 0.0;
  });

  solver::compute_diffusion_term(velocity, viscosity, diffusion);

  const double u_error = active_l2_error(diffusion.x, [viscosity](const double x, const double y, double) {
    return 2.0 * viscosity * std::cos(x) * std::sin(y);
  });
  const double v_error = active_l2_error(diffusion.y, [viscosity](const double x, const double y, double) {
    return -2.0 * viscosity * std::sin(x) * std::cos(y);
  });

  require(u_error < 2.0e-2, "diffusion u error too large");
  require(v_error < 2.0e-2, "diffusion v error too large");
}

struct TaylorGreenStepErrors {
  double velocity_error;
  double energy_error;
};

TaylorGreenStepErrors run_taylor_green_step_case(const int resolution) {
  const double domain_length = 2.0 * pi();
  const double viscosity = 0.01;
  const double spacing = domain_length / static_cast<double>(resolution);
  const solver::Grid grid{resolution, resolution, 1, spacing, spacing, 1.0, 1};
  const double dt = 0.1 * spacing * spacing;
  const solver::AdvectionOptions options{};

  solver::VelocityField velocity{grid};
  solver::VelocityField advection{grid};
  solver::VelocityField diffusion{grid};
  solver::VelocityField pressure_gradient{grid};
  solver::VelocityField updated_velocity{grid};
  solver::PressureField pressure{grid};

  fill_storage(velocity.x, [](const double x, const double y, double) {
    return -std::cos(x) * std::sin(y);
  });
  fill_storage(velocity.y, [](const double x, const double y, double) {
    return std::sin(x) * std::cos(y);
  });
  fill_storage(velocity.z, [](double, double, double) {
    return 0.0;
  });
  fill_storage(pressure, [](const double x, const double y, double) {
    return -0.25 * (std::cos(2.0 * x) + std::cos(2.0 * y));
  });

  updated_velocity = velocity;

  solver::compute_advection_term(velocity, options, advection);
  solver::compute_diffusion_term(velocity, viscosity, diffusion);
  solver::operators::compute_gradient(pressure, pressure_gradient);

  axpy_velocity(updated_velocity, advection, -dt);
  axpy_velocity(updated_velocity, pressure_gradient, -dt);
  axpy_velocity(updated_velocity, diffusion, dt);

  const double decay_velocity = std::exp(-2.0 * viscosity * dt);
  const double exact_energy = 0.25 * std::exp(-4.0 * viscosity * dt);
  const double numerical_energy = kinetic_energy(updated_velocity);

  const double u_error = active_l2_error(updated_velocity.x, [decay_velocity](const double x,
                                                                              const double y,
                                                                              double) {
    return -std::cos(x) * std::sin(y) * decay_velocity;
  });
  const double v_error = active_l2_error(updated_velocity.y, [decay_velocity](const double x,
                                                                              const double y,
                                                                              double) {
    return std::sin(x) * std::cos(y) * decay_velocity;
  });

  return TaylorGreenStepErrors{
      .velocity_error = std::sqrt(0.5 * (square(u_error) + square(v_error))),
      .energy_error = std::abs(numerical_energy - exact_energy),
  };
}

void test_taylor_green_step_behavior() {
  const TaylorGreenStepErrors coarse = run_taylor_green_step_case(16);
  const TaylorGreenStepErrors fine = run_taylor_green_step_case(32);

  require(coarse.velocity_error > fine.velocity_error,
          "Taylor-Green velocity error did not decrease on refinement");
  require(coarse.energy_error > fine.energy_error,
          "Taylor-Green energy error did not decrease on refinement");
  require(observed_order(coarse.velocity_error, fine.velocity_error) >= 1.5,
          "Taylor-Green velocity error did not show expected convergence");
  require(observed_order(coarse.energy_error, fine.energy_error) >= 1.5,
          "Taylor-Green energy error did not show expected convergence");
}

void test_bounded_advection_regression_case() {
  const int resolution_x = 64;
  const int resolution_y = 16;
  const double length = 1.0;
  const double dx = length / static_cast<double>(resolution_x);
  const double dy = length / static_cast<double>(resolution_y);
  const solver::Grid grid{resolution_x, resolution_y, 1, dx, dy, 1.0, 1};

  solver::VelocityField velocity{grid};
  solver::VelocityField advection{grid};
  solver::VelocityField updated{grid};
  solver::AdvectionOptions options{};
  const double dt = 0.4 * dx;

  fill_storage(velocity.x, [](double, double, double) {
    return 1.0;
  });
  fill_storage(velocity.y, [](double, double, double) {
    return 0.0;
  });
  fill_storage(velocity.z, [length](const double x, double, double) {
    const double wrapped_x = wrap_periodic(x, length);
    return wrapped_x < 0.5 ? 1.0 : 0.0;
  });

  updated = velocity;
  solver::compute_advection_term(velocity, options, advection);
  axpy_active(updated.z, advection.z, -dt);

  require(active_min(updated.z) >= -1.0e-12, "bounded TVD update undershot the minimum");
  require(active_max(updated.z) <= 1.0 + 1.0e-12, "bounded TVD update overshot the maximum");
}

void test_static_projection_preserves_quiescent_fluid() {
  const solver::Grid grid{16, 12, 1, 1.0 / 16.0, 1.0 / 12.0, 1.0, 1};
  const solver::BoundaryConditionSet boundary_conditions = solver::BoundaryConditionSet::cavity();
  const solver::ProjectionOptions options{
      .dt = 0.05,
      .density = 1.0,
      .poisson_max_iterations = 2000,
      .poisson_tolerance = 1.0e-12,
  };

  solver::VelocityField predicted{grid};
  solver::VelocityField corrected{grid};
  solver::PressureField pressure{grid};
  solver::ScalarField pressure_rhs{grid};
  predicted.fill(0.0);
  corrected.fill(0.0);
  pressure.fill(0.0);

  const solver::ProjectionDiagnostics diagnostics =
      solver::project_velocity(predicted, boundary_conditions, options, pressure, corrected, &pressure_rhs);

  require(diagnostics.pressure_solve.converged, "zero-RHS pressure solve should converge immediately");
  require(diagnostics.pressure_solve.iterations == 0, "quiescent projection should not iterate");
  require(diagnostics.rhs_l2 <= 1.0e-14, "quiescent projection RHS should remain zero");
  require(diagnostics.divergence_l2_before <= 1.0e-14,
          "quiescent predicted field should already be divergence free");
  require(diagnostics.divergence_l2_after <= 1.0e-14,
          "quiescent corrected field should remain divergence free");
  require(std::abs(diagnostics.pressure_mean) <= 1.0e-14, "quiescent pressure should remain zero");
  require(velocity_max_abs(corrected) <= 1.0e-12, "quiescent corrected velocity should remain zero");
}

void test_pure_neumann_projection_recovers_zero_mean_pressure() {
  const int resolution = 24;
  const solver::Grid grid{resolution, resolution, 1,
                          1.0 / static_cast<double>(resolution),
                          1.0 / static_cast<double>(resolution),
                          1.0,
                          1};
  const solver::BoundaryConditionSet boundary_conditions =
      solver::BoundaryConditionSet::all(solver::PhysicalBoundaryType::symmetry);
  const solver::ProjectionOptions options{
      .dt = 0.025,
      .density = 1.0,
      .poisson_max_iterations = 4000,
      .poisson_tolerance = 1.0e-12,
  };
  const solver::PressureBoundarySet pressure_boundaries =
      solver::derive_pressure_boundary_conditions(boundary_conditions);

  solver::PressureField expected_pressure{grid};
  fill_storage(expected_pressure, [](const double x, const double y, double) {
    const double sx = std::sin(pi() * x);
    const double sy = std::sin(pi() * y);
    return sx * sx * sy * sy;
  });
  solver::apply_pressure_boundary_conditions(pressure_boundaries, expected_pressure);

  solver::VelocityField predicted{grid};
  solver::operators::compute_gradient(expected_pressure, predicted);
  scale_velocity_active(predicted, options.dt);
  solver::apply_velocity_boundary_conditions(boundary_conditions, predicted);

  solver::PressureField pressure{grid};
  solver::VelocityField corrected{grid};
  pressure.fill(0.0);
  corrected.fill(0.0);

  const solver::ProjectionDiagnostics diagnostics =
      solver::project_velocity(predicted, boundary_conditions, options, pressure, corrected);

  const double expected_mean = active_mean_value(expected_pressure);
  const solver::IndexRange3D active = expected_pressure.layout().active_range();
  for(int k = active.k_begin; k < active.k_end; ++k) {
    for(int j = active.j_begin; j < active.j_end; ++j) {
      for(int i = active.i_begin; i < active.i_end; ++i) {
        expected_pressure(i, j, k) -= expected_mean;
      }
    }
  }
  solver::apply_pressure_boundary_conditions(pressure_boundaries, expected_pressure);

  const double pressure_error = active_l2_difference(pressure, expected_pressure);
  require(diagnostics.pressure_solve.converged, "pure-Neumann pressure solve should converge");
  require(diagnostics.pressure_solve.zero_mean_enforced,
          "pure-Neumann pressure solve should enforce zero mean");
  require(std::abs(diagnostics.pressure_mean) <= 1.0e-10,
          "pure-Neumann solve should remove the pressure null space");
  require(pressure_error <= 1.0e-8, "pressure recovery error is too large");
  require(velocity_max_abs(corrected) <= 1.0e-8, "projection should remove the gradient field");
  require(diagnostics.divergence_l2_after <= 1.0e-10,
          "projection should leave a divergence-free corrected field");
}

void test_lid_driven_cavity_config_loader_and_bc_subset() {
  const solver::LidDrivenCavityConfig config = solver::load_lid_driven_cavity_config(
      source_path("benchmarks/lid_driven_cavity_smoke.cfg"));
  require(config.nx == 32 && config.ny == 32, "smoke cavity config should load the grid size");
  require(std::abs(config.reynolds - 100.0) < 1.0e-12, "smoke cavity Reynolds number mismatch");
  require(!config.validate_reference, "smoke config should skip reference validation");

  const std::string summary = solver::describe(config);
  require(summary.find("validate_reference=false") != std::string::npos,
          "config description should report validation mode");

  const solver::BoundaryConditionSet boundary_conditions =
      solver::make_lid_driven_cavity_boundary_conditions(config);
  require(boundary_conditions[solver::BoundaryFace::x_min].type ==
              solver::PhysicalBoundaryType::no_slip_wall,
          "cavity x_min should remain a no-slip wall");
  require(boundary_conditions[solver::BoundaryFace::x_max].type ==
              solver::PhysicalBoundaryType::no_slip_wall,
          "cavity x_max should remain a no-slip wall");
  require(boundary_conditions[solver::BoundaryFace::y_min].type ==
              solver::PhysicalBoundaryType::no_slip_wall,
          "cavity y_min should remain a no-slip wall");
  require(boundary_conditions[solver::BoundaryFace::y_max].type ==
              solver::PhysicalBoundaryType::prescribed_velocity,
          "cavity lid should be a prescribed velocity boundary");
  require(std::abs(boundary_conditions[solver::BoundaryFace::y_max].velocity[0] -
                   config.lid_velocity) <= 1.0e-12,
          "cavity lid speed mismatch");
  require(boundary_conditions[solver::BoundaryFace::z_min].type ==
              solver::PhysicalBoundaryType::symmetry,
          "2D cavity z_min should map to symmetry");
  require(boundary_conditions[solver::BoundaryFace::z_max].type ==
              solver::PhysicalBoundaryType::symmetry,
          "2D cavity z_max should map to symmetry");
}

void test_lid_driven_cavity_smoke_run() {
  const solver::LidDrivenCavityConfig config = solver::load_lid_driven_cavity_config(
      source_path("benchmarks/lid_driven_cavity_smoke.cfg"));
  const solver::LidDrivenCavityResult result = solver::run_lid_driven_cavity(config);

  require(result.final_step.step == config.max_steps,
          "smoke cavity run should consume the fixed short-step budget");
  require(result.final_step.dt > 0.0, "smoke cavity timestep should be positive");
  require(result.final_step.max_cfl <= config.cfl_limit + 1.0e-12,
          "smoke cavity CFL should respect the configured ceiling");
  require(result.final_step.divergence_l2 <= 1.0e-10,
          "smoke cavity run should remain divergence controlled");
  require(result.final_step.pressure_iterations >= 0,
          "smoke cavity run should report pressure iterations");
  require(result.validation.reference_dataset.empty(),
          "smoke cavity run should skip benchmark validation");
  require(result.extrema.u_vertical_max > 0.0, "smoke cavity lid should induce positive u motion");
  require(result.extrema.u_vertical_min < 0.0,
          "smoke cavity recirculation should induce negative u motion");
}

void test_lid_driven_cavity_reference_validation_gate() {
  const solver::LidDrivenCavityReference reference = solver::ghia_re100_reference();
  solver::LidDrivenCavityResult passing{};
  passing.extrema.u_vertical_max = 0.0;
  passing.extrema.u_vertical_min = 0.0;
  passing.extrema.v_horizontal_max = 0.0;
  passing.extrema.v_horizontal_min = 0.0;
  passing.u_vertical_centerline.coordinate = {
      0.4000, reference.u_vertical_min_y, 0.5000, reference.u_vertical_max_y, 0.9900};
  passing.u_vertical_centerline.value = {
      -0.1800, reference.u_vertical_min, -0.20581, reference.u_vertical_max, 0.9000};
  passing.v_horizontal_centerline.coordinate = {
      0.2000, reference.v_horizontal_max_x, 0.5000, reference.v_horizontal_min_x, 0.9000};
  passing.v_horizontal_centerline.value = {
      0.1700, reference.v_horizontal_max, 0.05454, reference.v_horizontal_min, -0.1800};
  passing.final_step.divergence_l2 = 5.0e-11;

  const solver::LidDrivenCavityValidation validation =
      solver::validate_lid_driven_cavity_re100(passing);
  require(validation.pass, "reference-matching cavity result should pass validation");
  require(validation.reference_dataset == reference.dataset, "wrong reference dataset label");
  require(std::abs(validation.u_vertical_max_sample - reference.u_vertical_max) <= 1.0e-12,
          "validation should sample the vertical-max reference point");
  require(std::abs(validation.v_horizontal_min_sample - reference.v_horizontal_min) <= 1.0e-12,
          "validation should sample the horizontal-min reference point");

  solver::LidDrivenCavityResult failing = passing;
  failing.v_horizontal_centerline.value[1] *= 0.97;
  const solver::LidDrivenCavityValidation failed_validation =
      solver::validate_lid_driven_cavity_re100(failing);
  require(!failed_validation.pass, "3 percent extrema drift should fail the 2 percent gate");
}

}  // namespace

int main() {
  try {
    test_build_profile_is_locked();
    test_runtime_platform_is_supported();
    test_banner_contains_profile();
    test_grid_coordinates();
    test_pressure_layout_indexing();
    test_ghost_cell_access_and_boundary_ranges();
    test_memory_layout_and_alignment();
    test_cell_and_face_placement();
    test_double_precision_storage_contract();
    test_advection_options_and_cfl_diagnostic();
    test_projection_boundary_mapping();
    test_predictor_adi_preserves_quiescent_state();
    test_poisson_mgpcg_discrete_dirichlet_recovery();
    test_poisson_mgpcg_pure_neumann_zero_mean_recovery();
    test_manufactured_solution_convergence();
    test_diffusion_term_matches_scaled_laplacian();
    test_taylor_green_step_behavior();
    test_bounded_advection_regression_case();
    test_static_projection_preserves_quiescent_fluid();
    test_pure_neumann_projection_recovers_zero_mean_pressure();
    test_lid_driven_cavity_config_loader_and_bc_subset();
    test_lid_driven_cavity_smoke_run();
    test_lid_driven_cavity_reference_validation_gate();
  } catch(const std::exception& exception) {
    std::cerr << "solver_tests failed: " << exception.what() << '\n';
    return 1;
  }

  std::cout << "solver_tests passed" << '\n';
  return 0;
}
