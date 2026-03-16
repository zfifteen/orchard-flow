add_executable(solver_example
  example_main.cpp
)

target_link_libraries(solver_example
  PRIVATE
    solver_momentum
)

add_executable(solver_cavity
  lid_driven_cavity_main.cpp
)

target_link_libraries(solver_cavity
  PRIVATE
    solver_io
)

add_executable(solver_channel
  channel_flow_main.cpp
)

target_link_libraries(solver_channel
  PRIVATE
    solver_momentum
)

add_executable(solver_operator_mms
  operator_mms_main.cpp
)

target_link_libraries(solver_operator_mms
  PRIVATE
    solver_momentum
)

add_executable(solver_taylor_green
  taylor_green_main.cpp
)

target_link_libraries(solver_taylor_green
  PRIVATE
    solver_io
)

add_executable(solver_advection_profile
  advection_profile_main.cpp
)

target_link_libraries(solver_advection_profile
  PRIVATE
    solver_momentum
)

add_executable(solver_pressure_profile
  pressure_profile_main.cpp
)

target_link_libraries(solver_pressure_profile
  PRIVATE
    solver_momentum
)
