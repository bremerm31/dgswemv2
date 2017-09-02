#include "num_configuration.hpp"

#include "preprocessor/ADCIRC_reader/adcirc_format.hpp"
#include "preprocessor/mesh_metadata.hpp"

std::vector<std::vector<MeshMetaData>> partition(const MeshMetaData& mesh_meta,
                                                 const int num_partitions,
                                                 const int num_nodes,
                                                 const NumaConfiguration& numa_config);

bool check_partition(const MeshMetaData& mesh,
                     std::vector<std::vector<MeshMetaData>>& submeshes) {
  bool error_found;
  //2 checks are performed in this unit test
  //  1.) Check that all the elements can be found in one of the submeshes
  {
    std::unordered_set<uint> elts;
    for ( auto& e : mesh ) {
      elts.insert(e.first);
    }

    for ( auto& loc : submeshes ) {
      for ( auto& sm : loc ) {
        for ( auto& e : sm ) {
          elts.erase(e.first);
        }
      }
    }

    if ( !elts.empty() ) {
      error_found = true;
      std::cerr << "Error: Not all elements have been assigned to a submesh\n"
                << "       Specifically elements: {";
      for ( auto& e : elts ) {
        std::cerr << " " << e;
      }
      std::cerr << " } have not been assigned\n";
    }
  }


  //  2.) Check that all the edges between submeshes are set correctly

  return error_found;
}

//This unit test partitions a mesh and then stitches it back together again
int main(int argc, char** argv) {

  bool error_found {false};

  NumaConfiguration numa_config("default");

  AdcircFormat mesh1(argv[1]);
  MeshMetaData meshA(mesh1);

  std::cout << "Checking for 1 locality...\n";
  {
    std::vector<std::vector<MeshMetaData>> submeshes = partition(meshA, 4, 1, numa_config);
    error_found = error_found || check_partition(meshA, submeshes);
  }
  std::cout << "...done checking for 1 locality\n\n";
  std::cout << "Checking for multiple localities...\n";
  {
    std::vector<std::vector<MeshMetaData>> submeshes = partition(meshA, 4, 2, numa_config);
    error_found = error_found || check_partition(meshA, submeshes);
  }
  std::cout << "...done checking for multiple localities\n";

  if ( error_found ) {
    return 1;
  }

  return 0;
}
