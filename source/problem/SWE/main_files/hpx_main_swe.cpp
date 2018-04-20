#include <hpx/hpx.hpp>
#include <hpx/hpx_init.hpp>
#include <hpx/hpx_main.hpp>
#include <hpx/include/iostreams.hpp>
#include <hpx/include/components.hpp>
#include <hpx/include/lcos.hpp>

#include "../../../general_definitions.hpp"
#include "../swe_definitions.hpp"

#include "../function_files/swe_initial_condition_functions.hpp"
#include "../function_files/swe_source_functions.hpp"
#include "../function_files/swe_true_solution_functions.hpp"

#include "../swe_problem.hpp"
#include "../kernels_preprocessor/swe_kernels_preprocessor.hpp"
#include "../kernels_processor/swe_kernels_processor.hpp"
#include "../kernels_postprocessor/swe_kernels_postprocessor.hpp"

#include "../../../simulation/hpx_simulation.hpp"

DGSWEMV2_REGISTER_COMPONENTS(SWE::Problem);

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage\n"
                  << "    /path/to/DG_HYPER_SWE input_file\n";
        return 1;
    } else {
        return hpx::init(argc, argv);
    }
}

int hpx_main(int argc, char* argv[]) {
    std::string input_string = std::string(argv[1]);

    const std::vector<hpx::naming::id_type> localities = hpx::find_all_localities();

    std::vector<HPXSimulationClient<SWE::Problem>> simulation_clients;
    simulation_clients.reserve(localities.size());

    auto t1 = std::chrono::high_resolution_clock::now();
    for (hpx::naming::id_type const& locality : localities) {
        simulation_clients.emplace_back(
            hpx::new_<HPXSimulation<SWE::Problem>>(locality, input_string)
            );
    }

    std::vector<hpx::future<void>> run_futures;
    run_futures.reserve(simulation_clients.size());

    for (auto& sim_client : simulation_clients) {
        run_futures.push_back(sim_client.Run());
    }

    hpx::wait_all(run_futures);
    auto t2 = std::chrono::high_resolution_clock::now();

    std::cout << "Time Elapsed (in us): " << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count()
              << std::endl;

    return hpx::finalize();
}
