# Default-Policy Hotspots

| Symbol | Samples | Share | Category |
| --- | --- | --- | --- |
| solver::FieldLayout::is_storage_index(int, int, int) const | 133 | 0.347 | field_layout |
| solver::linsolve::(anonymous namespace)::neighbor_value(solver::PressureField const&, solver::PressureBoundarySet const&, int, int, int, solver::Axis, bool) | 28 | 0.073 | pressure_solve |
| solver::linsolve::solve_pressure_poisson(solver::ScalarField const&, solver::PressureBoundarySet const&, solver::ProjectionOptions const&, solver::PressureField&) | 16 | 0.042 | pressure_solve |
| solver::(anonymous namespace)::apply_face_field_boundary(solver::FaceField&, solver::BoundaryFace, solver::BoundaryCondition const&) | 14 | 0.037 | other |
| solver::(anonymous namespace)::solve_tridiagonal(std::__1::vector<double, std::__1::allocator<double>> const&, std::__1::vector<double, std::__1::allocator<double>> const&, std::__1::vector<double, std::__1::allocator<double>> const&, std::__1::vector<double, std::__1::allocator<double>> const&, std::__1::vector<double, std::__1::allocator<double>>&) | 11 | 0.029 | predictor_adi |
| solver::(anonymous namespace)::reconstruct_face_x(solver::StructuredField const&, int, int, int, double, solver::AdvectionOptions const&) | 10 | 0.026 | advection |
| solver::FieldLayout::storage_index_from_active(int, int, int) const | 10 | 0.026 | field_layout |
| solver::FieldLayout::index(int, int, int) const | 8 | 0.021 | field_layout |
| solver::(anonymous namespace)::axpy_velocity(solver::VelocityField&, solver::VelocityField const&, double) | 8 | 0.021 | other |
| solver::linsolve::(anonymous namespace)::apply_poisson_operator(solver::PressureField const&, solver::PressureBoundarySet const&, solver::PressureField&) | 7 | 0.018 | pressure_solve |
| solver::linsolve::(anonymous namespace)::active_mean(solver::StructuredField const&) | 7 | 0.018 | other |
| solver::(anonymous namespace)::compute_u_advection_2d(solver::VelocityField const&, solver::AdvectionOptions const&, solver::FaceField&) | 6 | 0.016 | advection |
