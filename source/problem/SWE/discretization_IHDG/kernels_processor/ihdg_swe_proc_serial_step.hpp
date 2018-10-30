#ifndef IHDG_SWE_PROC_SERIAL_STEP_HPP
#define IHDG_SWE_PROC_SERIAL_STEP_HPP

#include "ihdg_swe_kernels_processor.hpp"
#include "ihdg_swe_proc_serial_sol_glob_prob.hpp"

namespace SWE {
namespace IHDG {
template <typename SerialSimType>
void Problem::step_serial(SerialSimType* sim) {
    for (uint stage = 0; stage < sim->stepper.GetNumStages(); ++stage) {
        if (sim->parser.ParsingInput()) {
            sim->parser.ParseInput(sim->stepper, sim->discretization.mesh);
        }

        Problem::stage_serial(sim->stepper, sim->discretization);

        ++(sim->stepper);
    }

    if (sim->writer.WritingOutput()) {
        sim->writer.WriteOutput(sim->stepper, sim->discretization.mesh);
    }
}

template <typename StepperType, typename ProblemType>
void Problem::stage_serial(const StepperType& stepper, HDGDiscretization<ProblemType>& discretization) {
    discretization.mesh.CallForEachElement([&stepper, &discretization](auto& elt) {
        const uint stage = stepper.GetStage();

        auto& state      = elt.data.state[stage + 1];
        auto& state_prev = elt.data.state[stage];
        auto& internal   = elt.data.internal;

        state.q = state_prev.q;

        internal.q_prev_at_gp = elt.ComputeUgp(state_prev.q);
    });

    uint iter = 0;
    while (true) {
        iter++;

        /* Local Step */
        discretization.mesh.CallForEachElement([&stepper](auto& elt) { Problem::local_volume_kernel(stepper, elt); });

        discretization.mesh.CallForEachElement([&stepper](auto& elt) { Problem::local_source_kernel(stepper, elt); });

        discretization.mesh.CallForEachInterface(
            [&stepper](auto& intface) { Problem::local_interface_kernel(stepper, intface); });

        discretization.mesh.CallForEachBoundary(
            [&stepper](auto& bound) { Problem::local_boundary_kernel(stepper, bound); });

        discretization.mesh_skeleton.CallForEachEdgeInterface(
            [&stepper](auto& edge_int) { Problem::local_edge_interface_kernel(stepper, edge_int); });

        discretization.mesh_skeleton.CallForEachEdgeBoundary(
            [&stepper](auto& edge_bound) { Problem::local_edge_boundary_kernel(stepper, edge_bound); });

        /* Local Step */

        /* Global Step */
        discretization.mesh_skeleton.CallForEachEdgeInterface(
            [&stepper](auto& edge_int) { Problem::global_edge_interface_kernel(stepper, edge_int); });

        discretization.mesh_skeleton.CallForEachEdgeBoundary(
            [&stepper](auto& edge_bound) { Problem::global_edge_boundary_kernel(stepper, edge_bound); });
        /* Global Step */

        bool converged = Problem::serial_solve_global_problem(stepper, discretization);

        if (converged) {
            break;
        }

        if (iter == 100) {
            break;
        }
    }

    discretization.mesh.CallForEachElement([&stepper](auto& elt) {
        uint n_stages = stepper.GetNumStages();

        auto& state = elt.data.state;

        std::swap(state[0].q, state[n_stages].q);
    });

    discretization.mesh.CallForEachElement([&stepper](auto& elt) {
        bool nan_found = SWE::scrutinize_solution(stepper, elt);

        if (nan_found)
            abort();
    });
}
}
}

#endif