#define HERMES_REPORT_INFO
#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;

// This test makes sure that subdomains work correctly.

// Uniform polynomial degree of mesh elements.
const int P_INIT = 2;
// Number of initial uniform mesh refinements.
const int INIT_REF_NUM = 0;
// Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
// SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.
Hermes::MatrixSolverType matrix_solver_type = Hermes::SOLVER_UMFPACK;

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh_current, mesh_heat;
  Hermes::vector<Mesh*> meshes(&mesh_current, &mesh_heat);
  MeshReaderH2DXML mloader;
  mloader.load("problem-traverse.xml", meshes);

  // Perform initial mesh refinements (optional).
  for(int i = 0; i < INIT_REF_NUM; i++)
    for(unsigned int meshes_i = 0; meshes_i < meshes.size(); meshes_i++)
      meshes[meshes_i]->refine_all_elements();

  mloader.save("problem2.xml", meshes);
  mloader.load("problem2.xml", meshes);

  // Initialize essential boundary conditions.
  DefaultEssentialBCConst<double> bc_essential_current(Hermes::vector<std::string>("8", "9", "6", "3"), 0.0);
  EssentialBCs<double> bcs_current(&bc_essential_current);

  DefaultEssentialBCConst<double> bc_essential_heat(Hermes::vector<std::string>("2", "5", "1", "6", "4", "0", "7", "3"), 0.0);
  EssentialBCs<double> bcs_heat(&bc_essential_heat);

//  DefaultEssentialBCConst<double> bc_essential_complement(Hermes::vector<std::string>("Bottom Right", "Top Right", "Top Left", "Horizontal Left", "Vertical Bottom"), 0.0);
//  EssentialBCs<double> bcs_complement(&bc_essential_complement);

  // Create H1 spaces with default shapeset.
  H1Space<double> space_current(&mesh_current, &bcs_current, P_INIT);
  int ndof_current = space_current.get_num_dofs();

  H1Space<double> space_heat(&mesh_heat, &bcs_heat, P_INIT);
  int ndof_heat = space_heat.get_num_dofs();

//  H1Space<double> space_complement(&mesh_complement, &bcs_complement, P_INIT);
//  int ndof_complement = space_complement.get_num_dofs();

  info("ndofs current: %d\n", ndof_current);
  info("ndofs heat: %d\n", ndof_heat);
//  if(ndof_whole_domain == 9)
//  {
//    info("Success!");
//    return TEST_SUCCESS;
//  }
//  else
//  {
//    info("Failure!");
//    return TEST_FAILURE;
//  }

  return 0;
}
