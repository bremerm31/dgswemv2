#ifndef HPX_SIMULATION_UNIT_HPP
#define HPX_SIMULATION_UNIT_HPP

#include "../general_definitions.hpp"

#include "../preprocessor/input_parameters.hpp"
#include "../preprocessor/initialize_mesh.hpp"
#include "../communication/hpx_communicator.hpp"

#include "writer.hpp"
#include "load_balancer/base_model.hpp"
#include "load_balancer/abstract_load_balancer_factory.hpp"

template <typename ProblemType>
class HPXSimulationUnit
    : public  hpx::components::migration_support<
                  hpx::components::component_base<HPXSimulationUnit<ProblemType>>
              > {
  private:
/*    using BaseType = hpx::components::migration_support<
                         hpx::components::component_base<HPXSimulationUnit<ProblemType>>
                         >;*/

    uint state;
    Stepper stepper;

    Writer<ProblemType> writer;
    typename ProblemType::ProblemParserType parser;
    typename ProblemType::ProblemMeshType mesh;
    typename ProblemType::ProblemInputType problem_input;
    HPXCommunicator communicator;
    std::unique_ptr<LoadBalancer::SubmeshModel> submesh_model = nullptr;

  public:
    HPXSimulationUnit() = default;
    HPXSimulationUnit(const std::string& input_string, const uint locality_id, const uint submesh_id);
    HPXSimulationUnit(HPXSimulationUnit&& rhs) = default;

    HPXSimulationUnit& operator=(HPXSimulationUnit&& rhs)=default;

    hpx::future<void> Preprocessor();
    HPX_DEFINE_COMPONENT_ACTION(HPXSimulationUnit, Preprocessor, PreprocessorAction);

    void Launch();
    HPX_DEFINE_COMPONENT_ACTION(HPXSimulationUnit, Launch, LaunchAction);

    hpx::future<void> Step();
    HPX_DEFINE_COMPONENT_ACTION(HPXSimulationUnit, Step, StepAction);

    double ResidualL2();
    HPX_DEFINE_COMPONENT_ACTION(HPXSimulationUnit, ResidualL2, ResidualL2Action);


    //This function is part of me debugging. Ignore this for now. It doesn't get called in the code (at the moment)
    void SerializeAndUnserialize();
    HPX_DEFINE_COMPONENT_ACTION(HPXSimulationUnit, SerializeAndUnserialize, SerializeAndUnserializeAction);


    template <typename Archive>
    void save(Archive& ar, unsigned) const;

    template <typename Archive>
    void load(Archive& ar, unsigned);
    HPX_SERIALIZATION_SPLIT_MEMBER();

  private:
    void Parse();
    HPX_DEFINE_COMPONENT_ACTION(HPXSimulationUnit, Parse, ParseAction);

    hpx::future<void> Stage();
    HPX_DEFINE_COMPONENT_ACTION(HPXSimulationUnit, Stage, StageAction);

    hpx::future<void> Postprocessor();
    HPX_DEFINE_COMPONENT_ACTION(HPXSimulationUnit, Postprocessor, PostprocessorAction);

    void SwapStates();
    HPX_DEFINE_COMPONENT_ACTION(HPXSimulationUnit, SwapStates, SwapStatesAction);

};

template <typename ProblemType>
HPXSimulationUnit<ProblemType>::HPXSimulationUnit(const std::string& input_string, const uint locality_id, const uint submesh_id) {
    static_assert(std::is_move_constructible<HPXSimulationUnit<ProblemType>>::value, "Error: simulation unit must be move constructible");
    this->state = 0;
    InputParameters<typename ProblemType::ProblemInputType> input(input_string, locality_id, submesh_id);
    this->stepper = Stepper(input.rk.nstages, input.rk.order, input.dt, input.T_end);
    this->problem_input = input.problem_input;
    this->writer = Writer<ProblemType>(input, locality_id, submesh_id);
    this->parser = typename ProblemType::ProblemParserType(input, locality_id, submesh_id);
    this->mesh = typename ProblemType::ProblemMeshType(input.polynomial_order);
    this->communicator = HPXCommunicator(input.mesh_file_name.substr(0, input.mesh_file_name.find_last_of('.')) + ".dbmd",
                                         locality_id,
                                         submesh_id);
    this->submesh_model = LoadBalancer::AbstractFactory::create_submesh_model<ProblemType>(locality_id,submesh_id);
    //assert(this->submesh_model);
    if ( locality_id == 0 && submesh_id == 0 ) {
        std::cout << "Making submesh model on locality = 0 submesh_id =  0" << std::endl;
    }
    ProblemType::initialize_problem_parameters(this->problem_input);

    input.ReadMesh();

    this->mesh.SetMeshName(input.mesh_data.mesh_name);

    if (this->writer.WritingLog()) {
        this->writer.StartLog();

        this->writer.GetLogFile() << "Starting simulation with p=" << input.polynomial_order << " for "
                                  << input.mesh_data.mesh_name << " mesh" << std::endl << std::endl;
    }

    initialize_mesh<ProblemType, HPXCommunicator>(
        this->mesh, input.mesh_data, this->communicator, this->problem_input, this->writer);

        ProblemType::initialize_data_parallel_pre_send_kernel(this->mesh, input.mesh_data, this->problem_input);
}

template <typename ProblemType>
hpx::future<void> HPXSimulationUnit<ProblemType>::Preprocessor() {
    hpx::future<void> receive_future = this->communicator.ReceivePreprocAll(this->stepper.get_timestamp());

    this->communicator.SendPreprocAll(this->stepper.get_timestamp());

    return receive_future.then([this](auto&&) {
            ProblemType::initialize_data_parallel_post_receive_kernel(
            this->mesh, this->problem_input);
        });
}

template <typename ProblemType>
void HPXSimulationUnit<ProblemType>::Launch() {
    if (this->writer.WritingLog()) {
        this->writer.GetLogFile() << std::endl << "Launching Simulation!" << std::endl << std::endl;
    }

    if (this->writer.WritingOutput()) {
        this->writer.WriteFirstStep(this->stepper, this->mesh);
    }

    uint n_stages = this->stepper.get_num_stages();

    auto resize_data_container = [n_stages](auto& elt) { elt.data.resize(n_stages + 1); };

    this->mesh.CallForEachElement(resize_data_container);
}

template <typename ProblemType>
void HPXSimulationUnit<ProblemType>::Parse() {
    if (this->parser.ParsingInput()) {

        if (this->writer.WritingVerboseLog()) {
            this->writer.GetLogFile() << "Parsing input"
                                      << std::endl;
        }

        this->parser.ParseInput(this->stepper, this->mesh);
    }
}

template <typename ProblemType>
hpx::future<void> HPXSimulationUnit<ProblemType>::Stage() {
    if (this->writer.WritingVerboseLog()) {
        this->writer.GetLogFile() << "Current (time, stage): (" << this->stepper.get_t_at_curr_stage() << ','
                                  << this->stepper.get_stage() << ')' << std::endl << "Starting work before receive"
                                  << std::endl;
    }

    //this->Parse();

    if (this->writer.WritingVerboseLog()) {
        this->writer.GetLogFile() << "Finished Parse()" << std::endl;
    }

    auto distributed_boundary_send_kernel = [this](auto& dbound) {
        ProblemType::distributed_boundary_send_kernel(this->stepper, dbound);
    };

    auto volume_kernel = [this](auto& elt) { ProblemType::volume_kernel(this->stepper, elt); };

    auto source_kernel = [this](auto& elt) { ProblemType::source_kernel(this->stepper, elt); };

    auto interface_kernel = [this](auto& intface) { ProblemType::interface_kernel(this->stepper, intface); };

    auto boundary_kernel = [this](auto& bound) { ProblemType::boundary_kernel(this->stepper, bound); };

    if (this->writer.WritingVerboseLog()) {
        this->writer.GetLogFile() << "Exchanging data" << std::endl;
    }

    //hpx::future<void> receive_future = this->communicator.ReceiveAll(this->stepper.get_timestamp());

    //this->mesh.CallForEachDistributedBoundary(distributed_boundary_send_kernel);

    //this->communicator.SendAll(this->stepper.get_timestamp());

    if (this->writer.WritingVerboseLog()) {
        this->writer.GetLogFile() << "Starting work before receive" << std::endl;
    }

    auto assert_functor = [](auto& element)->bool {
        if ( element.master->phi_gp.size() > 0 ) { return true; }

        std::cout << "On locality " << hpx::find_here() << '\n'
        << "element " << element.GetID() << " master address: " << element.master << std::endl;
        return false;
    };
    this->mesh.CallForEachElement([assert_functor](auto& element) {
            assert(element.master);
            assert(assert_functor(element));
        });

    this->mesh.CallForEachElement(volume_kernel);

    this->mesh.CallForEachElement(source_kernel);

    this->mesh.CallForEachInterface(interface_kernel);

    this->mesh.CallForEachBoundary(boundary_kernel);

    if (this->writer.WritingVerboseLog()) {
        this->writer.GetLogFile() << "Finished work before receive" << std::endl
                                  << "Starting to wait on receive with timestamp: " << this->stepper.get_timestamp()
                                  << std::endl;
    }

//    return receive_future.then([this](auto&&) {
        if (this->writer.WritingVerboseLog()) {
            this->writer.GetLogFile() << "Starting work after receive" << std::endl;
        }

        auto distributed_boundary_kernel = [this](auto& dbound) {
            ProblemType::distributed_boundary_kernel(this->stepper, dbound);
        };

        auto update_kernel = [this](auto& elt) { ProblemType::update_kernel(this->stepper, elt); };

        auto scrutinize_solution_kernel = [this](auto& elt) {
            ProblemType::scrutinize_solution_kernel(this->stepper, elt);
        };

        //      this->mesh.CallForEachDistributedBoundary(distributed_boundary_kernel);

        this->mesh.CallForEachElement(update_kernel);

        this->mesh.CallForEachElement(scrutinize_solution_kernel);

        if (this->writer.WritingVerboseLog()) {
            this->writer.GetLogFile() << "Finished work after receive" << std::endl << std::endl;
        }
//        });
        return hpx::make_ready_future();
}

template <typename ProblemType>
hpx::future<void> HPXSimulationUnit<ProblemType>::Postprocessor() {
/*    if (this->writer.WritingVerboseLog()) {
        this->writer.GetLogFile() << "Exchanging postprocessor data" << std::endl;
    }

    hpx::future<void> receive_future = this->communicator.ReceivePostprocAll(this->stepper.get_timestamp());

    if (this->writer.WritingVerboseLog()) {
        this->writer.GetLogFile() << "postprocessor_parallel_pre_send_kernel" << std::endl;
    }

    ProblemType::postprocessor_parallel_pre_send_kernel(this->stepper, this->mesh);

    if (this->writer.WritingVerboseLog()) {
        this->writer.GetLogFile() << "sendpostprocall" << std::endl;
    }

    this->communicator.SendPostprocAll(this->stepper.get_timestamp());

    if (this->writer.WritingVerboseLog()) {
        this->writer.GetLogFile() << "Starting postprocessor work before receive" << std::endl;
    }

    ProblemType::postprocessor_parallel_pre_receive_kernel(this->stepper, this->mesh);

    if (this->writer.WritingVerboseLog()) {
        this->writer.GetLogFile() << "Finished postprocessor work before receive" << std::endl
                                  << "Starting to wait on postprocessor receive with timestamp: "
                                  << this->stepper.get_timestamp() << std::endl;
    }

    return receive_future.then([this](auto&&) {
        if (this->writer.WritingVerboseLog()) {
            this->writer.GetLogFile() << "Starting postprocessor work after receive" << std::endl;
        }

        ProblemType::postprocessor_parallel_post_receive_kernel(this->stepper, this->mesh);

        ++(this->stepper);

        if (this->writer.WritingVerboseLog()) {
            this->writer.GetLogFile() << "Finished postprocessor work after receive" << std::endl << std::endl;
        }
        });*/
    ++(this->stepper);
    return hpx::make_ready_future();
}

template <typename ProblemType>
void HPXSimulationUnit<ProblemType>::SwapStates() {
    auto swap_states_kernel = [this](auto& elt) { ProblemType::swap_states_kernel(this->stepper, elt); };

    this->mesh.CallForEachElement(swap_states_kernel);

    if (this->writer.WritingOutput()) {
        this->writer.WriteOutput(this->stepper, this->mesh);
    }
}

template <typename ProblemType>
hpx::future<void> HPXSimulationUnit<ProblemType>::Step() {
    auto assert_functor = [](auto& element)->bool {
        if ( element.master->phi_gp.size() > 0 ) { return true; }

        std::cout << "On locality " << hpx::find_here() << '\n'
                  << "element " << element.GetID() << " master address: " << element.master << std::endl;
        return false;
    };
    this->mesh.CallForEachElement([assert_functor](auto& element) {
            assert(element.master);
            assert(assert_functor(element)); //code throws here
        });


    hpx::future<void> step_future = hpx::make_ready_future();

    for (uint stage = 0; stage < this->stepper.get_num_stages(); stage++) {

        step_future = step_future.then([this](auto&&) {
                return this->Stage();
            });

        step_future = step_future.then([this](auto&&) {
                return this->Postprocessor();
                });
    }

    return step_future.then([this](auto&&) {
            assert(this->submesh_model);
            this->submesh_model->InStep(0,0);
            this->SwapStates();
            });
}

template <typename ProblemType>
double HPXSimulationUnit<ProblemType>::ResidualL2() {
    double residual_L2 = 0;

    auto compute_residual_L2_kernel = [this, &residual_L2](auto& elt) {
        residual_L2 += ProblemType::compute_residual_L2_kernel(this->stepper, elt);
    };

    this->mesh.CallForEachElement(compute_residual_L2_kernel);

    this->writer.GetLogFile() << "residual inner product: " << residual_L2 << std::endl;

    return residual_L2;
}

template <typename ProblemType>
void HPXSimulationUnit<ProblemType>::SerializeAndUnserialize() {
  std::cout << "serializing...\n";
  std::vector<char> buffer;
  hpx::serialization::output_archive oarchive(buffer);
  this->save(oarchive,0);

  stepper=Stepper();
  writer= std::move(Writer<ProblemType>());
  using ParserType  = typename ProblemType::ProblemParserType;
  parser=ParserType();
  using MeshType = typename ProblemType::ProblemMeshType;
  mesh = MeshType();
  using ProbInputType = typename ProblemType::ProblemInputType;
  problem_input = ProbInputType();
  std::unique_ptr<LoadBalancer::SubmeshModel> submesh_model = nullptr;

  hpx::serialization::input_archive iarchive(buffer);
  this->load(iarchive,0);
  std::cout << "...and unserializing" << std::endl;
}

template <typename ProblemType>
template <typename Archive>
void HPXSimulationUnit<ProblemType>::save(Archive& ar, unsigned) const {
    if (this->writer.WritingVerboseLog()) {
        std::pair<uint,uint> tag = submesh_model->get_tag();
        std::cout << tag.first << "," << tag.second << " departing from locality " << hpx::get_locality_id() << std::endl;
        this->writer.GetLogFile() << "Departing from locality " << hpx::get_locality_id() << std::endl;
    }

    ar & state & stepper & writer & parser & mesh & problem_input & submesh_model;
}

template <typename ProblemType>
template <typename Archive>
void HPXSimulationUnit<ProblemType>::load(Archive& ar, unsigned) {
    ar & state & stepper & writer & parser & mesh & problem_input & submesh_model;

    this->writer.StartLog();

    if (this->writer.WritingVerboseLog()) {
        std::pair<uint,uint> tag = submesh_model->get_tag();
        std::cout << tag.first << "," << tag.second << " arriving on locality " << hpx::get_locality_id() << std::endl;
        this->writer.GetLogFile() << "Arriving on locality " << hpx::get_locality_id() << std::endl;
    }

    this->mesh.CallForEachElement([](auto& element) {
            assert(element.master);
            assert(element.master->phi_gp.size() > 0);
            if ( element.GetID() == 330 ) {
                std::cout << "element " << element.GetID() << " master address: " << element.master << std::endl;
            }
        });


    initialize_mesh_interfaces_boundaries<ProblemType,HPXCommunicator>(mesh, communicator, writer);

    this->mesh.CallForEachElement([](auto& element) {
            assert(element.master);
            assert(element.master->phi_gp.size() > 0);
            if ( element.GetID() == 330 ) {
                std::cout << "element " << element.GetID() << " master address: " << element.master << std::endl;
            }

        });

}

template <typename ProblemType>
class HPXSimulationUnitClient
    : public hpx::components::client_base<HPXSimulationUnitClient<ProblemType>, HPXSimulationUnit<ProblemType>> {
  private:
    using BaseType = hpx::components::client_base<HPXSimulationUnitClient<ProblemType>, HPXSimulationUnit<ProblemType>>;

  public:
    HPXSimulationUnitClient()=default;
    HPXSimulationUnitClient(hpx::future<hpx::id_type>&& id) : BaseType(std::move(id)) {}
    HPXSimulationUnitClient(hpx::id_type&& id) : BaseType(std::move(id)) {}

    static constexpr const char* GetBasename() { return "Simulation_Unit_Client_"; }

    hpx::future<void> Preprocessor() {
        using ActionType = typename HPXSimulationUnit<ProblemType>::PreprocessorAction;
        return hpx::async<ActionType>(this->get_id());
    }

    hpx::future<void> Launch() {
        using ActionType = typename HPXSimulationUnit<ProblemType>::LaunchAction;
        return hpx::async<ActionType>(this->get_id());
    }

    hpx::future<void> Step() {
        using ActionType = typename HPXSimulationUnit<ProblemType>::StepAction;
        return hpx::async<ActionType>(this->get_id());
    }

    hpx::future<double> ResidualL2() {
        using ActionType = typename HPXSimulationUnit<ProblemType>::ResidualL2Action;
        return hpx::async<ActionType>(this->get_id());
    }

    hpx::future<void> SerializeAndUnserialize() {
        using ActionType = typename HPXSimulationUnit<ProblemType>::SerializeAndUnserializeAction;
        return hpx::async<ActionType>(this->get_id());
    }
};

#endif