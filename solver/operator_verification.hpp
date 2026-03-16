#pragma once

namespace solver {

struct OperatorManufacturedSolutionResult {
  int resolution = 0;
  double domain_length = 0.0;
  double gradient_error = 0.0;
  double divergence_error = 0.0;
  double laplacian_error = 0.0;
};

[[nodiscard]] OperatorManufacturedSolutionResult run_operator_manufactured_solution_case(
    int resolution);

}  // namespace solver
