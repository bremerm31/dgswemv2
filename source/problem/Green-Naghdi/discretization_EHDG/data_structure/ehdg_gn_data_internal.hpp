#ifndef EHDG_GN_DATA_INTERNAL_HPP
#define EHDG_GN_DATA_INTERNAL_HPP

#include "general_definitions.hpp"

namespace GN {
namespace EHDG {
struct Internal {
    Internal() = default;
    Internal(const uint ngp)
        : q_at_gp(GN::n_variables, ngp),
          aux_at_gp(GN::n_auxiliaries, ngp),
          Fx_at_gp(GN::n_variables, ngp),
          Fy_at_gp(GN::n_variables, ngp),
          source_at_gp(GN::n_variables, ngp),
          dbath_at_gp(GN::n_dimensions, ngp),
          tau_s_at_gp(GN::n_dimensions, ngp),
          dp_atm_at_gp(GN::n_dimensions, ngp),
          dtide_pot_at_gp(GN::n_dimensions, ngp) {}

    HybMatrix<double, GN::n_variables> q_at_gp;
    HybMatrix<double, GN::n_auxiliaries> aux_at_gp;

    HybMatrix<double, GN::n_variables> Fx_at_gp;
    HybMatrix<double, GN::n_variables> Fy_at_gp;

    HybMatrix<double, GN::n_variables> source_at_gp;
    HybMatrix<double, GN::n_dimensions> dbath_at_gp;
    HybMatrix<double, GN::n_dimensions> tau_s_at_gp;
    HybMatrix<double, GN::n_dimensions> dp_atm_at_gp;
    HybMatrix<double, GN::n_dimensions> dtide_pot_at_gp;
};
}
}

#endif