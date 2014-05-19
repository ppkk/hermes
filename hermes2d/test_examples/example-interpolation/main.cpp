#include "definitions.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Hermes2D::WeakFormsH1;

// This example shows how to solve a simple PDE that describes stationary
// heat transfer in an object consisting of two materials (aluminum and
// copper). The object is heated by constant volumetric heat sources
// (generated, for example, by a DC electric current). The temperature
// on the boundary is fixed. We will learn how to:
//
//   - load the mesh,
//   - perform initial refinements,
//   - create a H1 space over the mesh,
//   - define weak formulation,
//   - initialize matrix solver,
//   - assemble and solve the matrix system,
//   - output the solution and element orders in VTK format
//     (to be visualized, e.g., using Paraview),
//   - visualize the solution using Hermes' native OpenGL-based functionality.
//
// PDE: Poisson equation -div(LAMBDA grad u) - VOLUME_HEAT_SRC = 0.
//
// Boundary conditions: Dirichlet u(x, y) = FIXED_BDY_TEMP on the boundary.
//
// Geometry: L-Shape domain (see file domain.mesh).
//
// The following parameters can be changed:

const bool HERMES_VISUALIZATION = true;   // Set to "false" to suppress Hermes OpenGL visualization.
const bool VTK_VISUALIZATION = true;     // Set to "true" to enable VTK output.
const int P_INIT = 3;                     // Uniform polynomial degree of mesh elements.
const int INIT_REF_NUM = 3;               // Number of initial uniform mesh refinements.

// Problem parameters.
const double LAMBDA_AL = 236.0;            // Thermal cond. of Al for temperatures around 20 deg Celsius.
const double LAMBDA_CU = 386.0;            // Thermal cond. of Cu for temperatures around 20 deg Celsius.
const double VOLUME_HEAT_SRC = 5;        // Volume heat sources generated (for example) by electric current.
const double FIXED_BDY_TEMP = 20;        // Fixed temperature on the boundary.

int main(int argc, char* argv[])
{
  // Load the mesh.
  MeshSharedPtr mesh(new Mesh);
  Hermes::Hermes2D::MeshReaderH2D mloader;
  mloader.load("domain.mesh", mesh);
  mesh->refine_element_id(2);

  ScalarView s;
  BaseView<double> b;

  // Initialize essential boundary conditions.
  Hermes::Hermes2D::DefaultEssentialBCConst<double> bc_essential("Right", FIXED_BDY_TEMP);
  Hermes::Hermes2D::EssentialBCs<double> bcs(&bc_essential);

  // Initialize space.
  SpaceSharedPtr<double> space(new Hermes::Hermes2D::H1Space<double>(mesh, &bcs, 1));
  H1Space<double>* rawspace = dynamic_cast<H1Space<double>*>(space.get());
  b.show(space);

  // Initialize the weak formulation.
  DefaultWeakFormPoissonLinear<double> wf(HERMES_ANY, new Hermes::Hermes2DFunction<double>(VOLUME_HEAT_SRC));

  // Initialize the solution.
  MeshFunctionSharedPtr<double> sln(new Solution<double>);

  // Initialize linear solver.
  Hermes::Hermes2D::LinearSolver<double> linear_solver(&wf, space);

  // Solve.
  linear_solver.solve();
  Solution<double>::vector_to_solution(linear_solver.get_sln_vector(), space, sln);
  s.show(sln);

  //
  // Reference stuff.
  //
  ScalarView rs;
  BaseView<double> rb;

  // Construct globally refined reference mesh and setup reference space.
  Mesh::ReferenceMeshCreator ref_mesh_creator(mesh);
  MeshSharedPtr rmesh = ref_mesh_creator.create_ref_mesh();
  Space<double>::ReferenceSpaceCreator ref_space_creator(space, rmesh, 0);
  SpaceSharedPtr<double> rspace = ref_space_creator.create_ref_space();
  rb.show(rspace);

  // Initialize the solution.
  MeshFunctionSharedPtr<double> rsln(new Solution<double>);
  // Initialize linear solver.
  Hermes::Hermes2D::LinearSolver<double> rlinear_solver(&wf, rspace);
  rlinear_solver.solve();
  Solution<double>::vector_to_solution(rlinear_solver.get_sln_vector(), rspace, rsln);
  rs.show(rsln);
  View::wait();

  //
  // Interpolated.
  //
  double* interpolated_vector;
  MeshFunctionSharedPtr<double> isln(new Solution<double>);
  VertexBasedInterpolation<double>::interpolate(rspace, rlinear_solver.get_sln_vector(), space, interpolated_vector);
  Solution<double>::vector_to_solution(interpolated_vector, space, isln);
  ScalarView is;
  is.show(isln);

  View::wait();
  return 0;
}