#include "problem.h"
#include "definitions.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;

const int P_INIT = 3;                     // Uniform polynomial degree of mesh elements.
const int INIT_REF_NUM = 1;               // Number of initial uniform mesh refinements.

// Problem parameters.


MeshFunctionSharedPtr<double> solve_problem(ProblemDefinition definition, MeshSharedPtr mesh)
{
    // Initialize essential boundary conditions.
    Hermes::Hermes2D::DefaultEssentialBCConst<double> bc_essential_ground(definition.bc_ground, 0);
    Hermes::Hermes2D::DefaultEssentialBCConst<double> bc_essential_potential(definition.bc_potential, definition.POTENTIAL);
    Hermes::Hermes2D::EssentialBCs<double> bcs(Hermes::vector<EssentialBoundaryCondition<double> *> (&bc_essential_ground, &bc_essential_potential));

    SpaceSharedPtr<double> space(new Hermes::Hermes2D::H1Space<double>(mesh, &bcs, P_INIT));
    std::cout << "Ndofs: " << space->get_num_dofs() << std::endl;

    CustomWeakFormPoisson wf(definition);
    MeshFunctionSharedPtr<double> sln(new Solution<double>);
    Hermes::Hermes2D::LinearSolver<double> linear_solver(&wf, space);

    try
    {
        linear_solver.solve();
        double* sln_vector = linear_solver.get_sln_vector();

        Hermes::Hermes2D::Solution<double>::vector_to_solution(sln_vector, space, sln);
        return sln;
    }
    catch (Exceptions::Exception& e)
    {
        std::cout << e.info();
    }
    catch (std::exception& e)
    {
        std::cout << e.what();
    }
}

MeshFunctionSharedPtr<double> solve_permitivity(ProblemDefinition definition, MeshSharedPtr mesh)
{
    SpaceSharedPtr<double> space_perm(new Hermes::Hermes2D::H1Space<double>(mesh, nullptr, 3));
    CustomWeakFormPermitivity wf_perm(definition);
    MeshFunctionSharedPtr<double> sln_perm(new Solution<double>);
    Hermes::Hermes2D::LinearSolver<double> linear_solver_perm(&wf_perm, space_perm);

    try
    {
        linear_solver_perm.solve();
        double* sln_vector_perm = linear_solver_perm.get_sln_vector();
        Hermes::Hermes2D::Solution<double>::vector_to_solution(sln_vector_perm, space_perm, sln_perm);
        return sln_perm;
    }
    catch (Exceptions::Exception& e)
    {
        std::cout << e.info();
    }
    catch (std::exception& e)
    {
        std::cout << e.what();
    }
}

double calc_integral(MeshFunctionSharedPtr<double> sln, ProblemDefinition definition)
{
    MyVolumetricIntegralCalculator integralCalculator(sln, definition);
    double* result = new double[1];
    result = integralCalculator.calculate(HERMES_ANY);
    double res = result[0];
    delete[] result;
    return res;
}

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

    MeshFunctionSharedPtr<double> sln = solve_problem(definition, mesh);
    MeshFunctionSharedPtr<double> sln_perm = solve_permitivity(definition, mesh);

    std::cout << "Integral is " << calc_integral(sln, definition) << std::endl;

    // Visualize the solution.
    Hermes::Hermes2D::Views::ScalarView viewS("Solution", new Hermes::Hermes2D::Views::WinGeom(0, 0, 1500, 700));
    Hermes::Hermes2D::Views::ScalarView viewP("Permitivity", new Hermes::Hermes2D::Views::WinGeom(0, 700, 1500, 700));
    //    Hermes::Hermes2D::Views::OrderView viewSp("Space", new Hermes::Hermes2D::Views::WinGeom(0, 600, 1200, 600));
    //    viewSp.show(space);
    viewS.show(sln);
    viewP.show(sln_perm);
    viewS.wait_for_close();
    return 0;
}
