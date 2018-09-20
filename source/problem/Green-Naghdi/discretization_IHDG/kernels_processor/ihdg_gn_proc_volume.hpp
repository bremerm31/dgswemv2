#ifndef IHDG_GN_PROC_VOLUME_HPP
#define IHDG_GN_PROC_VOLUME_HPP

namespace GN {
namespace IHDG {
template <typename ElementType>
void Problem::local_swe_volume_kernel(const RKStepper& stepper, ElementType& elt) {
    const uint stage = stepper.GetStage();

    auto& state    = elt.data.state[stage];
    auto& internal = elt.data.internal;

    internal.q_at_gp = elt.ComputeUgp(state.q);

    row(internal.aux_at_gp, GN::Auxiliaries::h) =
        row(internal.q_at_gp, GN::Variables::ze) + row(internal.aux_at_gp, GN::Auxiliaries::bath);

    auto u = vec_cw_div(row(internal.q_at_gp, GN::Variables::qx), row(internal.aux_at_gp, GN::Auxiliaries::h));
    auto v = vec_cw_div(row(internal.q_at_gp, GN::Variables::qy), row(internal.aux_at_gp, GN::Auxiliaries::h));

    auto uuh = vec_cw_mult(u, row(internal.q_at_gp, GN::Variables::qx));
    auto vvh = vec_cw_mult(v, row(internal.q_at_gp, GN::Variables::qy));
    auto uvh = vec_cw_mult(u, row(internal.q_at_gp, GN::Variables::qy));
    auto pe  = Global::g *
              (0.5 * vec_cw_mult(row(internal.q_at_gp, GN::Variables::ze), row(internal.q_at_gp, GN::Variables::ze)) +
               vec_cw_mult(row(internal.q_at_gp, GN::Variables::ze), row(internal.aux_at_gp, GN::Auxiliaries::bath)));

    row(internal.Fx_at_gp, GN::Variables::ze) = row(internal.q_at_gp, GN::Variables::qx);
    row(internal.Fx_at_gp, GN::Variables::qx) = uuh + pe;
    row(internal.Fx_at_gp, GN::Variables::qy) = uvh;

    row(internal.Fy_at_gp, GN::Variables::ze) = row(internal.q_at_gp, GN::Variables::qy);
    row(internal.Fy_at_gp, GN::Variables::qx) = uvh;
    row(internal.Fy_at_gp, GN::Variables::qy) = vvh + pe;

    state.rhs =
        elt.IntegrationDPhi(GlobalCoord::x, internal.Fx_at_gp) + elt.IntegrationDPhi(GlobalCoord::y, internal.Fy_at_gp);
}

template <typename ElementType>
void Problem::local_dc_volume_kernel(const RKStepper& stepper, ElementType& elt) {
    auto& internal = elt.data.internal;

    // at this point h_at_gp
    // has been calculated in derivatives kernel

    // set kernels up
    double bx, by;

    for (uint gp = 0; gp < elt.data.get_ngp_internal(); ++gp) {
        bx = internal.dbath_at_gp(GlobalCoord::x, gp);
        by = internal.dbath_at_gp(GlobalCoord::y, gp);

        column(internal.w1_w1_kernel_at_gp, gp) = IdentityVector<double>(GN::n_dimensions);

        internal.w1_w1_kernel_at_gp(GN::n_dimensions * GlobalCoord::x + GlobalCoord::x, gp) +=
            NDParameters::alpha * bx * bx;
        internal.w1_w1_kernel_at_gp(GN::n_dimensions * GlobalCoord::x + GlobalCoord::y, gp) +=
            NDParameters::alpha * bx * by;
        internal.w1_w1_kernel_at_gp(GN::n_dimensions * GlobalCoord::y + GlobalCoord::x, gp) +=
            NDParameters::alpha * bx * by;
        internal.w1_w1_kernel_at_gp(GN::n_dimensions * GlobalCoord::y + GlobalCoord::y, gp) +=
            NDParameters::alpha * by * by;
    }

    set_constant(internal.w1_w2_kernel_at_gp, -NDParameters::alpha / 3.0);

    row(internal.w2_w1_kernel_at_gp, GlobalCoord::x) = power(row(internal.aux_at_gp, GN::Auxiliaries::h), -1.0);
    row(internal.w2_w1_kernel_at_gp, GlobalCoord::y) = power(row(internal.aux_at_gp, GN::Auxiliaries::h), -1.0);

    internal.w2_w2_kernel_at_gp = power(row(internal.aux_at_gp, GN::Auxiliaries::h), -3.0);

    for (uint dof_i = 0; dof_i < elt.data.get_ndof(); ++dof_i) {
        for (uint dof_j = 0; dof_j < elt.data.get_ndof(); ++dof_j) {
            submatrix(internal.w1_w1,
                      GN::n_dimensions * dof_i,
                      GN::n_dimensions * dof_j,
                      GN::n_dimensions,
                      GN::n_dimensions) =
                reshape<double, GN::n_dimensions>(elt.IntegrationPhiPhi(dof_i, dof_j, internal.w1_w1_kernel_at_gp));

            internal.w1_w2(GN::n_dimensions * dof_i + GlobalCoord::x, dof_j) =
                elt.IntegrationPhiDPhi(dof_i, GlobalCoord::x, dof_j, row(internal.w1_w2_kernel_at_gp, GlobalCoord::x));
            internal.w1_w2(GN::n_dimensions * dof_i + GlobalCoord::y, dof_j) =
                elt.IntegrationPhiDPhi(dof_i, GlobalCoord::y, dof_j, row(internal.w1_w2_kernel_at_gp, GlobalCoord::y));

            internal.w2_w1(dof_i, GN::n_dimensions * dof_j + GlobalCoord::x) =
                elt.IntegrationPhiDPhi(dof_j, GlobalCoord::x, dof_i, row(internal.w2_w1_kernel_at_gp, GlobalCoord::x));
            internal.w2_w1(dof_i, GN::n_dimensions * dof_j + GlobalCoord::y) =
                elt.IntegrationPhiDPhi(dof_j, GlobalCoord::y, dof_i, row(internal.w2_w1_kernel_at_gp, GlobalCoord::y));

            internal.w2_w2(dof_i, dof_j) = elt.IntegrationPhiPhi(dof_i, dof_j, internal.w2_w2_kernel_at_gp);
        }
    }

    set_constant(internal.w1_w1_kernel_at_gp, 0.0);

    row(internal.w1_w1_kernel_at_gp, GN::n_dimensions * GlobalCoord::x + GlobalCoord::x) =
        -NDParameters::alpha / 2.0 *
        vec_cw_mult(row(internal.aux_at_gp, GN::Auxiliaries::h), row(internal.dbath_at_gp, GlobalCoord::x));
    row(internal.w1_w1_kernel_at_gp, GN::n_dimensions * GlobalCoord::x + GlobalCoord::y) =
        -NDParameters::alpha / 2.0 *
        vec_cw_mult(row(internal.aux_at_gp, GN::Auxiliaries::h), row(internal.dbath_at_gp, GlobalCoord::y));
    row(internal.w1_w1_kernel_at_gp, GN::n_dimensions * GlobalCoord::y + GlobalCoord::x) =
        -NDParameters::alpha / 2.0 *
        vec_cw_mult(row(internal.aux_at_gp, GN::Auxiliaries::h), row(internal.dbath_at_gp, GlobalCoord::x));
    row(internal.w1_w1_kernel_at_gp, GN::n_dimensions * GlobalCoord::y + GlobalCoord::y) =
        -NDParameters::alpha / 2.0 *
        vec_cw_mult(row(internal.aux_at_gp, GN::Auxiliaries::h), row(internal.dbath_at_gp, GlobalCoord::y));

    set_constant(internal.w1_w2_kernel_at_gp, 0.0);

    row(internal.w1_w2_kernel_at_gp, GlobalCoord::x) =
        -NDParameters::alpha / 2.0 *
        vec_cw_div(row(internal.dbath_at_gp, GlobalCoord::x), row(internal.aux_at_gp, GN::Auxiliaries::h));
    row(internal.w1_w2_kernel_at_gp, GlobalCoord::y) =
        -NDParameters::alpha / 2.0 *
        vec_cw_div(row(internal.dbath_at_gp, GlobalCoord::y), row(internal.aux_at_gp, GN::Auxiliaries::h));

    for (uint dof_i = 0; dof_i < elt.data.get_ndof(); ++dof_i) {
        for (uint dof_j = 0; dof_j < elt.data.get_ndof(); ++dof_j) {
            internal.w1_w1(GN::n_dimensions * dof_i + GlobalCoord::x, GN::n_dimensions * dof_j + GlobalCoord::x) +=
                elt.IntegrationPhiDPhi(
                    dof_j,
                    GlobalCoord::x,
                    dof_i,
                    row(internal.w1_w1_kernel_at_gp, GN::n_dimensions * GlobalCoord::x + GlobalCoord::x));

            internal.w1_w1(GN::n_dimensions * dof_i + GlobalCoord::x, GN::n_dimensions * dof_j + GlobalCoord::y) +=
                elt.IntegrationPhiDPhi(
                    dof_j,
                    GlobalCoord::x,
                    dof_i,
                    row(internal.w1_w1_kernel_at_gp, GN::n_dimensions * GlobalCoord::x + GlobalCoord::y));

            internal.w1_w1(GN::n_dimensions * dof_i + GlobalCoord::y, GN::n_dimensions * dof_j + GlobalCoord::x) +=
                elt.IntegrationPhiDPhi(
                    dof_j,
                    GlobalCoord::y,
                    dof_i,
                    row(internal.w1_w1_kernel_at_gp, GN::n_dimensions * GlobalCoord::y + GlobalCoord::x));

            internal.w1_w1(GN::n_dimensions * dof_i + GlobalCoord::y, GN::n_dimensions * dof_j + GlobalCoord::y) +=
                elt.IntegrationPhiDPhi(
                    dof_j,
                    GlobalCoord::y,
                    dof_i,
                    row(internal.w1_w1_kernel_at_gp, GN::n_dimensions * GlobalCoord::y + GlobalCoord::y));

            internal.w1_w2(GN::n_dimensions * dof_i + GlobalCoord::x, dof_j) +=
                elt.IntegrationPhiPhi(dof_i, dof_j, row(internal.w1_w2_kernel_at_gp, GlobalCoord::x));
            internal.w1_w2(GN::n_dimensions * dof_i + GlobalCoord::y, dof_j) +=
                elt.IntegrationPhiPhi(dof_i, dof_j, row(internal.w1_w2_kernel_at_gp, GlobalCoord::y));
        }
    }
}
}
}

#endif
