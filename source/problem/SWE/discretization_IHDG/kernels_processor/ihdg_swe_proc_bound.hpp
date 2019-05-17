#ifndef IHDG_SWE_PROC_BOUND_HPP
#define IHDG_SWE_PROC_BOUND_HPP

namespace SWE {
namespace IHDG {
template <typename BoundaryType>
void Problem::init_boundary_kernel(const ProblemStepperType& stepper, BoundaryType& bound) {
    if (stepper.GetOrder() == 2) {
        const uint stage = stepper.GetStage();

        auto& state_prev = bound.data.state[stage];
        auto& boundary   = bound.data.boundary[bound.bound_id];

        for ( uint var = 0; var < SWE::n_variables; ++var ) {
            boundary.q_at_gp[var] = bound.ComputeUgp(state_prev.q[var]);
        }
    }
}

template <typename BoundaryType>
void Problem::local_boundary_kernel(const ProblemStepperType& stepper, BoundaryType& bound) {
    const uint stage = stepper.GetStage();

    auto& state    = bound.data.state[stage + 1];
    auto& boundary = bound.data.boundary[bound.bound_id];

    for ( uint var = 0; var < SWE::n_variables; ++var ) {
        boundary.q_at_gp[var] = bound.ComputeUgp(state.q[var]);
    }
}
}
}

#endif
