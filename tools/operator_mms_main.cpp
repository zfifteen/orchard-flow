#include "solver/operator_verification.hpp"

#include <iostream>
#include <stdexcept>
#include <string>

int main(int argc, char** argv) {
  try {
    const int resolution = argc > 1 ? std::stoi(argv[1]) : 32;
    const solver::OperatorManufacturedSolutionResult result =
        solver::run_operator_manufactured_solution_case(resolution);

    std::cout << "resolution: " << result.resolution << '\n';
    std::cout << "domain_length: " << result.domain_length << '\n';
    std::cout << "gradient_l2_error: " << result.gradient_error << '\n';
    std::cout << "divergence_l2_error: " << result.divergence_error << '\n';
    std::cout << "laplacian_l2_error: " << result.laplacian_error << '\n';
    return 0;
  } catch(const std::exception& exception) {
    std::cerr << "solver_operator_mms failed: " << exception.what() << '\n';
    return 1;
  }
}
