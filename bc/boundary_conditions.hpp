#pragma once

#include "core/fields.hpp"

#include <array>
#include <string>

namespace solver {

enum class PhysicalBoundaryType : int {
  no_slip_wall = 0,
  prescribed_velocity = 1,
  symmetry = 2,
  fixed_pressure = 3,
  periodic = 4,
};

enum class PressureBoundaryType : int {
  neumann = 0,
  dirichlet = 1,
  periodic = 2,
};

struct BoundaryCondition {
  PhysicalBoundaryType type = PhysicalBoundaryType::no_slip_wall;
  std::array<double, 3> velocity{0.0, 0.0, 0.0};
  double pressure = 0.0;
  double pressure_gradient = 0.0;
};

struct BoundaryConditionSet {
  std::array<BoundaryCondition, 6> faces{};

  [[nodiscard]] static BoundaryConditionSet all(PhysicalBoundaryType type);
  [[nodiscard]] static BoundaryConditionSet cavity();

  [[nodiscard]] BoundaryCondition& operator[](BoundaryFace face) noexcept {
    return faces[static_cast<std::size_t>(face)];
  }

  [[nodiscard]] const BoundaryCondition& operator[](BoundaryFace face) const noexcept {
    return faces[static_cast<std::size_t>(face)];
  }
};

struct PressureBoundaryCondition {
  PressureBoundaryType type = PressureBoundaryType::neumann;
  double value = 0.0;
  double gradient = 0.0;
};

struct PressureBoundarySet {
  std::array<PressureBoundaryCondition, 6> faces{};

  [[nodiscard]] PressureBoundaryCondition& operator[](BoundaryFace face) noexcept {
    return faces[static_cast<std::size_t>(face)];
  }

  [[nodiscard]] const PressureBoundaryCondition& operator[](BoundaryFace face) const noexcept {
    return faces[static_cast<std::size_t>(face)];
  }
};

std::string to_string(PhysicalBoundaryType type);
std::string to_string(PressureBoundaryType type);

PressureBoundarySet derive_pressure_correction_boundary_conditions(
    const BoundaryConditionSet& boundary_conditions);

void apply_velocity_boundary_conditions(const BoundaryConditionSet& boundary_conditions,
                                        VelocityField& velocity);

void apply_pressure_boundary_conditions(const PressureBoundarySet& boundary_conditions,
                                        PressureField& pressure);

// Fills total-pressure ghosts for the physical pressure field. Velocity-like
// normal gradient source terms are sampled from the matching wall-normal
// component for Neumann-style physical boundaries.
void apply_total_pressure_boundary_conditions(const BoundaryConditionSet& boundary_conditions,
                                              const VelocityField& normal_pressure_gradient_source,
                                              PressureField& pressure_total);

}  // namespace solver
