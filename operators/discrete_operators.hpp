#pragma once

#include "core/fields.hpp"

namespace solver::operators {

void compute_gradient(const PressureField& pressure, VelocityField& gradient);
void compute_divergence(const VelocityField& velocity, ScalarField& divergence);
void compute_laplacian(const StructuredField& input, StructuredField& output);
void compute_laplacian(const VelocityField& input, VelocityField& output);

}  // namespace solver::operators

