#include "problem.h"
#include "definitions.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;



const int P_INIT = 3;                     // Uniform polynomial degree of mesh elements.
const int INIT_REF_NUM = 1;               // Number of initial uniform mesh refinements.

// Problem parameters.

int main(int argc, char* argv[])
{
    Hermes::vector<int> profile_coarse(0,0,0,0,0,0,0,0,0,5);
    int active_electrode = 7;
    double eps_rel_material = 5;

    ProblemConfiguration configuration(profile_coarse, active_electrode);
    ProblemDefinition_1 definition(eps_rel_material, configuration);

    // Load the mesh.
    MeshSharedPtr mesh(new Mesh);
    Hermes::Hermes2D::MeshReaderH2DXML mloader;
    mloader.load(definition.mesh_name, mesh);

    // Refine all elements, do it INIT_REF_NUM-times.
    for (unsigned int i = 0; i < INIT_REF_NUM; i++)
        mesh->refine_all_elements();

    // Initialize essential boundary conditions.
    Hermes::Hermes2D::DefaultEssentialBCConst<double> bc_essential_ground(definition.bc_ground, 0);
    Hermes::Hermes2D::DefaultEssentialBCConst<double> bc_essential_potential(definition.bc_potential, definition.POTENTIAL);
    Hermes::Hermes2D::EssentialBCs<double> bcs(Hermes::vector<EssentialBoundaryCondition<double> *> (&bc_essential_ground, &bc_essential_potential));

    // Initialize space.
    SpaceSharedPtr<double> space(new Hermes::Hermes2D::H1Space<double>(mesh, &bcs, P_INIT));
    SpaceSharedPtr<double> space_perm(new Hermes::Hermes2D::H1Space<double>(mesh, nullptr, 3));

    std::cout << "Ndofs: " << space->get_num_dofs() << std::endl;


    // Initialize the weak formulation.
    CustomWeakFormPoisson wf(definition);

    CustomWeakFormPermitivity wf_perm(definition);

    // Initialize the solution.
    MeshFunctionSharedPtr<double> sln(new Solution<double>);
    MeshFunctionSharedPtr<double> sln_perm(new Solution<double>);

    // Initialize linear solver.
    Hermes::Hermes2D::LinearSolver<double> linear_solver(&wf, space);
    Hermes::Hermes2D::LinearSolver<double> linear_solver_perm(&wf_perm, space_perm);

    // Solve the linear problem.
    try
    {
        linear_solver.solve();
        linear_solver_perm.solve();

        // Get the solution vector.
        double* sln_vector = linear_solver.get_sln_vector();
        double* sln_vector_perm = linear_solver_perm.get_sln_vector();

        // Translate the solution vector into the previously initialized Solution.
        Hermes::Hermes2D::Solution<double>::vector_to_solution(sln_vector, space, sln);
        Hermes::Hermes2D::Solution<double>::vector_to_solution(sln_vector_perm, space_perm, sln_perm);

        MyVolumetricIntegralCalculator integralCalculator(sln, definition);
        double* result = new double[1];
        result = integralCalculator.calculate(HERMES_ANY);
        std::cout << "Integral is " << result[0] << std::endl;

        // Visualize the solution.
        Hermes::Hermes2D::Views::ScalarView viewS("Solution", new Hermes::Hermes2D::Views::WinGeom(0, 0, 1500, 700));
        Hermes::Hermes2D::Views::ScalarView viewP("Permitivity", new Hermes::Hermes2D::Views::WinGeom(0, 700, 1500, 700));
        Hermes::Hermes2D::Views::OrderView viewSp("Space", new Hermes::Hermes2D::Views::WinGeom(0, 600, 1200, 600));
        viewSp.show(space);
        viewS.show(sln);
        viewP.show(sln_perm);
        viewS.wait_for_close();
    }
    catch (Exceptions::Exception& e)
    {
        std::cout << e.info();
    }
    catch (std::exception& e)
    {
        std::cout << e.what();
    }
    return 0;
}
