#ifndef GN_PRE_CREATE_EDGE_BOUND_HPP
#define GN_PRE_CREATE_EDGE_BOUND_HPP

namespace GN {
template <typename ProblemType>
void create_edge_boundaries(typename ProblemType::ProblemMeshType& mesh,
                            typename ProblemType::ProblemMeshSkeletonType& mesh_skeleton,
                            typename ProblemType::ProblemWriterType& writer) {
    using BoundaryTypeLand = typename std::tuple_element<0, typename ProblemType::ProblemBoundaryTypes>::type;
    using BoundaryTypeTide = typename std::tuple_element<1, typename ProblemType::ProblemBoundaryTypes>::type;
    using BoundaryTypeFlow = typename std::tuple_element<2, typename ProblemType::ProblemBoundaryTypes>::type;

    using EdgeBoundaryTypeLand = typename std::tuple_element<0, typename ProblemType::ProblemEdgeBoundaryTypes>::type;
    using EdgeBoundaryTypeTide = typename std::tuple_element<1, typename ProblemType::ProblemEdgeBoundaryTypes>::type;
    using EdgeBoundaryTypeFlow = typename std::tuple_element<2, typename ProblemType::ProblemEdgeBoundaryTypes>::type;

    uint old_edge_land = mesh_skeleton.GetNumberEdgeBoundaries();

    mesh.template CallForEachBoundaryOfType<BoundaryTypeLand>(
        [&mesh_skeleton](auto& bound) { mesh_skeleton.template CreateEdgeBoundary<EdgeBoundaryTypeLand>(bound); });

    if (writer.WritingLog()) {
        writer.GetLogFile() << "Number of land edges: " << mesh_skeleton.GetNumberEdgeBoundaries() - old_edge_land
                            << std::endl;
    }

    uint old_edge_tide = mesh_skeleton.GetNumberEdgeBoundaries();

    mesh.template CallForEachBoundaryOfType<BoundaryTypeTide>(
        [&mesh_skeleton](auto& bound) { mesh_skeleton.template CreateEdgeBoundary<EdgeBoundaryTypeTide>(bound); });

    if (writer.WritingLog()) {
        writer.GetLogFile() << "Number of tide edges: " << mesh_skeleton.GetNumberEdgeBoundaries() - old_edge_tide
                            << std::endl;
    }

    uint old_edge_flow = mesh_skeleton.GetNumberEdgeBoundaries();

    mesh.template CallForEachBoundaryOfType<BoundaryTypeFlow>(
        [&mesh_skeleton](auto& bound) { mesh_skeleton.template CreateEdgeBoundary<EdgeBoundaryTypeFlow>(bound); });

    if (writer.WritingLog()) {
        writer.GetLogFile() << "Number of flow edges: " << mesh_skeleton.GetNumberEdgeBoundaries() - old_edge_flow
                            << std::endl;
    }
}
}

#endif
