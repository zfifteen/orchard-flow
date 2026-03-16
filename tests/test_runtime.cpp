#include "core/fields.hpp"
#include "core/grid.hpp"
#include "core/runtime.hpp"

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
  } catch(const std::exception& exception) {
    std::cerr << "solver_tests failed: " << exception.what() << '\n';
    return 1;
  }

  std::cout << "solver_tests passed" << '\n';
  return 0;
}
