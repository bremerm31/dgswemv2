#ifndef IHDG_SWE_PROC_EDGE_DBOUND_HPP
#define IHDG_SWE_PROC_EDGE_DBOUND_HPP

#include "problem/SWE/problem_flux/swe_flux.hpp"

namespace SWE {
namespace IHDG {
template <typename EdgeDistributedType>
void Problem::init_edge_distributed_kernel(const ProblemStepperType& stepper, EdgeDistributedType& edge_dbound) {
    if (stepper.GetOrder() == 2) {
        auto& edge_state    = edge_dbound.edge_data.edge_state;
        auto& edge_internal = edge_dbound.edge_data.edge_internal;

        auto& internal = edge_dbound.boundary.data.internal;
        auto& boundary = edge_dbound.boundary.data.boundary[edge_dbound.boundary.bound_id];

        auto& q_hat_at_gp    = edge_internal.q_hat_at_gp;
        auto& aux_hat_at_gp  = edge_internal.aux_hat_at_gp;
        auto& surface_normal = edge_dbound.boundary.surface_normal;

        edge_internal.q_hat_at_gp = edge_dbound.ComputeUgp(edge_state.q_hat);

        row(edge_internal.aux_hat_at_gp, SWE::Auxiliaries::h) =
            row(edge_internal.q_hat_at_gp, SWE::Variables::ze) +
            row(edge_internal.aux_hat_at_gp, SWE::Auxiliaries::bath);

        /* Compute fluxes at boundary states */

        SWE::get_Fn(q_hat_at_gp, aux_hat_at_gp, surface_normal, boundary.F_hat_at_gp);

        /* Add stabilization parameter terms */

        SWE::get_tau_LF(q_hat_at_gp, aux_hat_at_gp, surface_normal, edge_internal.tau);

        for (uint gp = 0; gp < edge_dbound.edge_data.get_ngp(); ++gp) {
            for ( uint var = 0; var < SWE::n_variables; ++var ) {
                for ( uint w = 0; w < SWE::n_variables; ++w ) {
                    boundary.F_hat_at_gp[var][gp] +=
                        edge_internal.tau[gp](var,w) * (boundary.q_at_gp[w][gp] - edge_internal.q_hat_at_gp(w, gp));
                }
            }
        }

        for (uint dof_i = 0; dof_i < edge_dbound.boundary.data.get_ndof(); ++dof_i) {
            for ( uint var = 0; var < SWE::n_variables; ++var ) {
                internal.rhs_prev[var + SWE::n_variables * dof_i] +=
                    -edge_dbound.boundary.IntegrationPhi(dof_i, boundary.F_hat_at_gp[var]);
            }
        }
    }
}

template <typename EdgeDistributedType>
void Problem::local_edge_distributed_kernel(const ProblemStepperType& stepper, EdgeDistributedType& edge_dbound) {
    auto& edge_state    = edge_dbound.edge_data.edge_state;
    auto& edge_internal = edge_dbound.edge_data.edge_internal;

    auto& internal = edge_dbound.boundary.data.internal;
    auto& boundary = edge_dbound.boundary.data.boundary[edge_dbound.boundary.bound_id];

    auto& q_hat_at_gp    = edge_internal.q_hat_at_gp;
    auto& aux_hat_at_gp  = edge_internal.aux_hat_at_gp;
    auto& surface_normal = edge_dbound.boundary.surface_normal;

    edge_internal.q_hat_at_gp = edge_dbound.ComputeUgp(edge_state.q_hat);

    row(edge_internal.aux_hat_at_gp, SWE::Auxiliaries::h) =
        row(edge_internal.q_hat_at_gp, SWE::Variables::ze) + row(edge_internal.aux_hat_at_gp, SWE::Auxiliaries::bath);

    /* Compute fluxes at boundary states */

    SWE::get_Fn(q_hat_at_gp, aux_hat_at_gp, surface_normal, boundary.F_hat_at_gp);

    SWE::get_dFn_dq(q_hat_at_gp, aux_hat_at_gp, surface_normal, boundary.dF_hat_dq_hat_at_gp);

    /* Add stabilization parameter terms */

    SWE::get_tau_LF(q_hat_at_gp, aux_hat_at_gp, surface_normal, edge_internal.tau);
    SWE::get_dtau_dze_LF(q_hat_at_gp, aux_hat_at_gp, surface_normal, edge_internal.dtau_dze);
    SWE::get_dtau_dqx_LF(q_hat_at_gp, aux_hat_at_gp, surface_normal, edge_internal.dtau_dqx);
    SWE::get_dtau_dqy_LF(q_hat_at_gp, aux_hat_at_gp, surface_normal, edge_internal.dtau_dqy);

    StatVector<double, SWE::n_variables> del_q;
    StatMatrix<double, SWE::n_variables, SWE::n_variables> dtau_delq;

    for (uint gp = 0; gp < edge_dbound.edge_data.get_ngp(); ++gp) {
        for ( uint var = 0; var < SWE::n_variables; ++var ) {
            del_q[var] = boundary.q_at_gp[var][gp] - edge_internal.q_hat_at_gp(var, gp);
        }

        column(dtau_delq, SWE::Variables::ze) = edge_internal.dtau_dze[gp] * del_q;
        column(dtau_delq, SWE::Variables::qx) = edge_internal.dtau_dqx[gp] * del_q;
        column(dtau_delq, SWE::Variables::qy) = edge_internal.dtau_dqy[gp] * del_q;

        for ( uint var = 0; var < SWE::n_variables; ++var ) {
            for ( uint w = 0; w < SWE::n_variables; ++w ) {
                boundary.F_hat_at_gp[var][gp] += edge_internal.tau[gp](var,w) * del_q[w];
            }
        }
        column(boundary.dF_hat_dq_at_gp, gp) = flatten<double>(edge_internal.tau[gp]);
        column(boundary.dF_hat_dq_hat_at_gp, gp) += flatten<double>(dtau_delq - edge_internal.tau[gp]);
    }

    double theta = stepper.GetTheta();

    for (uint dof_i = 0; dof_i < edge_dbound.boundary.data.get_ndof(); ++dof_i) {
        for (uint dof_j = 0; dof_j < edge_dbound.boundary.data.get_ndof(); ++dof_j) {
            submatrix(internal.delta_local,
                      SWE::n_variables * dof_i,
                      SWE::n_variables * dof_j,
                      SWE::n_variables,
                      SWE::n_variables) +=
                reshape<double, SWE::n_variables>(
                    (1.0 - theta) * edge_dbound.boundary.IntegrationPhiPhi(dof_j, dof_i, boundary.dF_hat_dq_at_gp));
        }

        for ( uint var = 0; var < SWE::n_variables; ++var ) {
            internal.rhs_local[ var + SWE::n_variables * dof_i] +=
                -(1.0 - theta) * edge_dbound.boundary.IntegrationPhi(dof_i, boundary.F_hat_at_gp[var]);
        }
    }

    for (uint dof_i = 0; dof_i < edge_dbound.boundary.data.get_ndof(); ++dof_i) {
        for (uint dof_j = 0; dof_j < edge_dbound.edge_data.get_ndof(); ++dof_j) {
            submatrix(boundary.delta_hat_local,
                      SWE::n_variables * dof_i,
                      SWE::n_variables * dof_j,
                      SWE::n_variables,
                      SWE::n_variables) =
                reshape<double, SWE::n_variables>(
                    (1.0 - theta) * edge_dbound.IntegrationPhiLambda(dof_i, dof_j, boundary.dF_hat_dq_hat_at_gp));
        }
    }
}

template <typename EdgeDistributedType>
void Problem::global_edge_distributed_kernel(const ProblemStepperType& stepper, EdgeDistributedType& edge_dbound) {
    auto& edge_internal = edge_dbound.edge_data.edge_internal;

    auto& boundary = edge_dbound.boundary.data.boundary[edge_dbound.boundary.bound_id];

    edge_dbound.boundary.boundary_condition.ComputeGlobalKernels(edge_dbound);

    for (uint dof_i = 0; dof_i < edge_dbound.edge_data.get_ndof(); ++dof_i) {
        for (uint dof_j = 0; dof_j < edge_dbound.boundary.data.get_ndof(); ++dof_j) {
            submatrix(boundary.delta_global,
                      SWE::n_variables * dof_i,
                      SWE::n_variables * dof_j,
                      SWE::n_variables,
                      SWE::n_variables) =
                reshape<double, SWE::n_variables>(
                    edge_dbound.IntegrationPhiLambda(dof_j, dof_i, boundary.delta_global_kernel_at_gp));
        }
    }

    for (uint dof_i = 0; dof_i < edge_dbound.edge_data.get_ndof(); ++dof_i) {
        for (uint dof_j = 0; dof_j < edge_dbound.edge_data.get_ndof(); ++dof_j) {
            submatrix(edge_internal.delta_hat_global,
                      SWE::n_variables * dof_i,
                      SWE::n_variables * dof_j,
                      SWE::n_variables,
                      SWE::n_variables) =
                reshape<double, SWE::n_variables>(
                    edge_dbound.IntegrationLambdaLambda(dof_j, dof_i, edge_internal.delta_hat_global_kernel_at_gp));
        }

        subvector(edge_internal.rhs_global, SWE::n_variables * dof_i, SWE::n_variables) =
            -edge_dbound.IntegrationLambda(dof_i, edge_internal.rhs_global_kernel_at_gp);
    }
}
}
}

#endif