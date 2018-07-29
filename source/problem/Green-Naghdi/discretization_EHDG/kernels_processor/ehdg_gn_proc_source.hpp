#ifndef EHDG_GN_PROC_SOURCE_HPP
#define EHDG_GN_PROC_SOURCE_HPP

namespace GN {
namespace EHDG {
template <typename ElementType>
void Problem::local_source_kernel(const RKStepper& stepper, ElementType& elt) {
    const uint stage = stepper.GetStage();

    auto& state    = elt.data.state[stage];
    auto& internal = elt.data.internal;
    // auto& source   = elt.data.source;

    double t = stepper.GetTimeAtCurrentStage();

    // note we assume that the values at gauss points have already been computed
    // compute contribution of hydrostatic pressure
    set_constant(row(internal.source_at_gp, GN::Variables::ze), 0.0);

    row(internal.source_at_gp, GN::Variables::qx) =
        Global::g *
        cwise_multiplication(row(internal.dbath_at_gp, GlobalCoord::x), row(internal.q_at_gp, GN::Variables::ze));

    row(internal.source_at_gp, GN::Variables::qy) =
        Global::g *
        cwise_multiplication(row(internal.dbath_at_gp, GlobalCoord::y), row(internal.q_at_gp, GN::Variables::ze));

    if (GN::SourceTerms::function_source) {
        auto source_u = [t](Point<2>& pt) { return GN::source_u(t, pt); };

        internal.source_at_gp += elt.ComputeFgp(source_u);
    }

    /*if (GN::SourceTerms::bottom_friction) {
        double Cf = GN::SourceTerms::Cf;

        for (uint gp = 0; gp < elt.data.get_ngp_internal(); ++gp) {
            // compute bottom friction contribution
            double u = internal.q_at_gp(GN::Variables::qx, gp) / internal.aux_at_gp(GN::Auxiliaries::h, gp);
            double v = internal.q_at_gp(GN::Variables::qy, gp) / internal.aux_at_gp(GN::Auxiliaries::h, gp);

            // compute manning friction factor
            if (source.manning) {
                Cf = source.g_manning_n_sq / std::pow(internal.aux_at_gp(GN::Auxiliaries::h, gp), 1.0 / 3.0);
                if (Cf < GN::SourceTerms::Cf)
                    Cf = GN::SourceTerms::Cf;
            }

            double bottom_friction_stress = Cf * std::hypot(u, v) / internal.aux_at_gp(GN::Auxiliaries::h, gp);

            internal.source_at_gp(GN::Variables::qx, gp) -=
                bottom_friction_stress * internal.q_at_gp(GN::Variables::qx, gp);
            internal.source_at_gp(GN::Variables::qy, gp) -=
                bottom_friction_stress * internal.q_at_gp(GN::Variables::qy, gp);
        }
    }

    if (GN::SourceTerms::meteo_forcing) {
        // elt.ComputeNodalUgp(source.tau_s, internal.tau_s_at_gp);
        // elt.ComputeNodalDUgp(source.p_atm, internal.dp_atm_at_gp);

        for (uint gp = 0; gp < elt.data.get_ngp_internal(); ++gp) {
            // compute surface friction contribution
            internal.source_at_gp(GN::Variables::qx, gp) +=
                internal.tau_s_at_gp(GlobalCoord::x, gp) / Global::rho_water;
            internal.source_at_gp(GN::Variables::qy, gp) +=
                internal.tau_s_at_gp(GlobalCoord::y, gp) / Global::rho_water;

            // compute atmospheric pressure contribution
            internal.source_at_gp(GN::Variables::qx, gp) -=
                internal.aux_at_gp(GN::Auxiliaries::sp, gp) * internal.aux_at_gp(GN::Auxiliaries::h, gp) *
                internal.dp_atm_at_gp(GlobalCoord::x, gp) / Global::rho_water;
            internal.source_at_gp(GN::Variables::qy, gp) -= internal.aux_at_gp(GN::Auxiliaries::h, gp) *
                                                             internal.dp_atm_at_gp(GlobalCoord::y, gp) /
                                                             Global::rho_water;
        }
    }

    if (GN::SourceTerms::tide_potential) {
        // elt.ComputeNodalDUgp(source.tide_pot, internal.dtide_pot_at_gp);

        for (uint gp = 0; gp < elt.data.get_ngp_internal(); ++gp) {
            // compute tide potential contribution
            internal.source_at_gp(GN::Variables::qx, gp) +=
                Global::g * internal.aux_at_gp(GN::Auxiliaries::sp, gp) *
                internal.aux_at_gp(GN::Auxiliaries::h, gp) * internal.dtide_pot_at_gp(GlobalCoord::x, gp);
            internal.source_at_gp(GN::Variables::qy, gp) += Global::g *
                                                             internal.aux_at_gp(GN::Auxiliaries::h, gp) *
                                                             internal.dtide_pot_at_gp(GlobalCoord::y, gp);
        }
    }

    if (GN::SourceTerms::coriolis) {
        for (uint gp = 0; gp < elt.data.get_ngp_internal(); ++gp) {
            // compute coriolis contribution
            internal.source_at_gp(GN::Variables::qx, gp) +=
                source.coriolis_f * internal.q_at_gp(GN::Variables::qy, gp);
            internal.source_at_gp(GN::Variables::qy, gp) -=
                source.coriolis_f * internal.q_at_gp(GN::Variables::qx, gp);
        }
    }*/

    state.rhs += elt.IntegrationPhi(internal.source_at_gp);
}
}
}

#endif