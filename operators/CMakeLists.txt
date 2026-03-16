add_library(solver_operators STATIC
  discrete_operators.cpp
)

target_include_directories(solver_operators
  PUBLIC
    ${PROJECT_SOURCE_DIR}
)

target_link_libraries(solver_operators
  PUBLIC
    solver_core
)

