#define HERMES_REPORT_INFO
#include "hermes2d.h"

#include "definitions.h"

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
const double NEWTON_TOL = 1e-4;
// Maximum allowed number of Newton iterations.
const int NEWTON_MAX_ITER = 100;

int main(int argc, char* argv[])
{
//  // Load the mesh.
//  Mesh mesh;
//  Hermes::vector<Mesh*> meshes;
//  meshes.push_back(&mesh);
//  MeshReaderH2DXML mloader;
//  mloader.load("bad_ori.xml", meshes);

//  // Perform initial mesh refinements (optional).
//  for(int i = 0; i < INIT_REF_NUM; i++)
//    for(unsigned int meshes_i = 0; meshes_i < meshes.size(); meshes_i++)
//      meshes[meshes_i]->refine_all_elements();

//  mloader.save("problem2.xml", meshes);
//  mloader.load("problem2.xml", meshes);

//  // Initialize essential boundary conditions.
//  DefaultEssentialBCConst<double> bc_essential(Hermes::vector<std::string>("20", "21", "22", "23"), 0.0);
//  EssentialBCs<double> bcs(&bc_essential);

//  // Create H1 spaces with default shapeset.
//  H1Space<double> space(&mesh, &bcs, P_INIT);
//  int ndof = space.get_num_dofs();

//  //info("ndofs current: %d\n", ndof);

//  // Initialize weak formulation.
//  CustomWeakFormCurrent wf;

//  // Initialize the FE problem.
//  DiscreteProblem<double> dp(&wf, &space);

//  // Initialize the Newton solver.
//  NewtonSolver<double> newton(&dp, matrix_solver_type);
//  double *coeff_vec = new double[ndof];

//    // Perform Newton's iteration.
//    info("Solving nonlinear problem:");
//    bool verbose = true;
//    // Perform Newton's iteration and translate the resulting coefficient vector into previous time level solutions.
//    newton.set_verbose_output(verbose);
//    try
//    {
//      newton.solve(coeff_vec, NEWTON_TOL, NEWTON_MAX_ITER);
//    }
//    catch(Hermes::Exceptions::Exception e)
//    {
//      e.printMsg();
//      error("Newton's iteration failed for current.");
//    };

//    Solution<double> solution;
//    Solution<double>::vector_to_solution(newton.get_sln_vector(), &space, &solution);



////    {
////      Hermes::vector<Solution<double> *> tmp(&xvel_prev_time, &yvel_prev_time, &p_prev_time, &temperature_prev_time);
////      Solution<double>::vector_to_solutions(newton.get_sln_vector(), Hermes::vector<const Space<double> *>(&xvel_space,
////          &yvel_space, &p_space, &temperature_space), tmp);
////    }

////    // Show the solution at the end of time step.
////    sprintf(title, "Velocity [m/s], time %g s", current_time);
////    vview.set_title(title);
////    vview.show(&xvel_prev_time, &yvel_prev_time);
////    sprintf(title, "Pressure [Pa], time %g s", current_time);
////    pview.set_title(title);
////    pview.show(&p_prev_time);
////    sprintf(title, "Temperature [C], time %g s", current_time);
////    tempview.set_title(title);
////    tempview.show(&temperature_prev_time, Views::HERMES_EPS_HIGH);
////  }

//  delete [] coeff_vec;

//  // Wait for all views to be closed.
//  Views::View::wait();
//  return 0;


////  if(ndof_whole_domain == 9)
////  {
////    info("Success!");
////    return TEST_SUCCESS;
////  }
////  else
////  {
////    info("Failure!");
////    return TEST_FAILURE;
////  }

//  return 0;
}
