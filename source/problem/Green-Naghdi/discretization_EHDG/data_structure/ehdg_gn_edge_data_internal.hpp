#ifndef EHDG_GN_EDGE_DATA_INTERNAL_HPP
#define EHDG_GN_EDGE_DATA_INTERNAL_HPP

#include "general_definitions.hpp"

namespace GN {
namespace EHDG {
struct EdgeInternal {
    EdgeInternal() = default;
    EdgeInternal(const uint ngp)
        : q_hat_at_gp(GN::n_variables, ngp),
          aux_hat_at_gp(GN::n_auxiliaries, ngp),
          q_init_at_gp(GN::n_variables, ngp),
          delta_hat_global_kernel_at_gp(GN::n_variables * GN::n_variables, ngp),
          rhs_global_kernel_at_gp(GN::n_variables, ngp) {}

    HybMatrix<double, GN::n_variables> q_hat_at_gp;
    HybMatrix<double, GN::n_auxiliaries> aux_hat_at_gp;
    HybMatrix<double, GN::n_variables> q_init_at_gp;

    HybMatrix<double, GN::n_variables * GN::n_variables> delta_hat_global_kernel_at_gp;
    HybMatrix<double, GN::n_variables> rhs_global_kernel_at_gp;

    DynMatrix<double> delta_hat_global;
    DynVector<double> rhs_global;
};
}
}

#endif