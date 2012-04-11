#define HERMES_REPORT_INFO
#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;

// This test makes sure that subdomains work correctly.

// Uniform polynomial degree of mesh elements.
const int P_INIT = 2;
// Number of initial uniform mesh refinements.
const int INIT_REF_NUM = 1;
// Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
// SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.
Hermes::MatrixSolverType matrix_solver_type = Hermes::SOLVER_UMFPACK;

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh_whole_domain, mesh_bottom_left_corner, mesh_complement;
  Hermes::vector<Mesh*> meshes(&mesh_whole_domain, &mesh_bottom_left_corner, &mesh_complement);
  MeshReaderH2DXML mloader;
  mloader.load("problem-multi.xml", meshes);

  // Perform initial mesh refinements (optional).
  for(int i = 0; i < INIT_REF_NUM; i++)
    for(unsigned int meshes_i = 0; meshes_i < meshes.size(); meshes_i++)
      meshes[meshes_i]->refine_all_elements();

  mloader.save("problem2.xml", meshes);
  mloader.load("problem2.xml", meshes);

  // Initialize essential boundary conditions.
  DefaultEssentialBCConst<double> bc_essential_whole_domain(Hermes::vector<std::string>("0", "1", "2", "3", "4", "5"), 0.0);
  EssentialBCs<double> bcs_whole_domain(&bc_essential_whole_domain);

//  DefaultEssentialBCConst<double> bc_essential_bottom_left_corner(Hermes::vector<std::string>("Bottom Left", "Horizontal Left"), 0.0);
//  EssentialBCs<double> bcs_bottom_left_corner(&bc_essential_bottom_left_corner);

//  DefaultEssentialBCConst<double> bc_essential_complement(Hermes::vector<std::string>("Bottom Right", "Top Right", "Top Left", "Horizontal Left", "Vertical Bottom"), 0.0);
//  EssentialBCs<double> bcs_complement(&bc_essential_complement);

  // Create H1 spaces with default shapeset.
  H1Space<double> space_whole_domain(&mesh_whole_domain, &bcs_whole_domain, P_INIT);
  int ndof_whole_domain = space_whole_domain.get_num_dofs();

//  H1Space<double> space_bottom_left_corner(&mesh_bottom_left_corner, &bcs_bottom_left_corner, P_INIT);
//  int ndof_bottom_left_corner = space_bottom_left_corner.get_num_dofs();

//  H1Space<double> space_complement(&mesh_complement, &bcs_complement, P_INIT);
//  int ndof_complement = space_complement.get_num_dofs();

  info("ndofs: %d\n", ndof_whole_domain);
  if(ndof_whole_domain == 9)
  {
    info("Success!");
    return TEST_SUCCESS;
  }
  else
  {
    info("Failure!");
    return TEST_FAILURE;
  }

  return 0;
}
