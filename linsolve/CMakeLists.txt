add_library(solver_linsolve STATIC
  poisson_solver.cpp
)

target_include_directories(solver_linsolve
  PUBLIC
    ${PROJECT_SOURCE_DIR}
)

target_link_libraries(solver_linsolve
  PUBLIC
    solver_core
)
