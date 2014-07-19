#include "problem.h"
#include "definitions.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;

const int P_INIT = 3;                     // Uniform polynomial degree of mesh elements.
const int INIT_REF_NUM = 1;               // Number of initial uniform mesh refinements.

// Problem parameters.


MeshFunctionSharedPtr<double> solve_problem(ProblemDefinition definition, Perms perms, MeshSharedPtr mesh)
{
    // Initialize essential boundary conditions.
    Hermes::Hermes2D::DefaultEssentialBCConst<double> bc_essential_ground(definition.bc_labels_ground, 0);
    Hermes::Hermes2D::DefaultEssentialBCConst<double> bc_essential_potential(definition.bc_labels_potential, definition.POTENTIAL);
    Hermes::Hermes2D::EssentialBCs<double> bcs(Hermes::vector<EssentialBoundaryCondition<double> *> (&bc_essential_ground, &bc_essential_potential));

    SpaceSharedPtr<double> space(new Hermes::Hermes2D::H1Space<double>(mesh, &bcs, P_INIT));
    std::cout << "Ndofs: " << space->get_num_dofs() << std::endl;

    CustomWeakFormPoisson wf(definition, perms);
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

MeshFunctionSharedPtr<double> solve_permitivity(ProblemDefinition definition, Perms perms, MeshSharedPtr mesh)
{
    SpaceSharedPtr<double> space_perm(new Hermes::Hermes2D::H1Space<double>(mesh, nullptr, 3));
    CustomWeakFormPermitivity wf_perm(definition, perms);
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

double calc_integral_energy(MeshFunctionSharedPtr<double> sln, ProblemDefinition definition, Perms perms)
{
    EnergyIntegralCalculator integralCalculator(sln, definition, perms);
    double* result = new double[1];
    result = integralCalculator.calculate(HERMES_ANY);
    double res = result[0];
    delete[] result;
    return res;
}

double calc_integral_grad_u_grad_v(MeshFunctionSharedPtr<double> sln1, MeshFunctionSharedPtr<double> sln2, Hermes::vector<std::string> area)
{
    GradUGradVIngegralCalculator integralCalculator(Hermes::vector<MeshFunctionSharedPtr<double> >(sln1, sln2));
    double* result = new double[1];
    result = integralCalculator.calculate(area);
    double res = result[0];
    delete[] result;
    return res;
}

double calc_integral_u_f(MeshFunctionSharedPtr<double> sln, double coeff)
{
    U_times_f_IngegralCalculator integralCalculator(sln, coeff);
    double* result = new double[1];
    result = integralCalculator.calculate(HERMES_ANY);
    double res = result[0];
    delete[] result;
    return res;
}

void update_fn_1d(PGDSolutions pgd_solutions, ProblemDefinition definition, Perms perms)
{
    // the last pair solutions/parameters is the actually calculated
    // the last parameter will be replaced by the now calculated (next iteration)
    assert(pgd_solutions.solutions.size() == pgd_solutions.parameters.size());
    double a = 0;
    double b = 0;   //a*eps + b
    double c = 0;   //---------
    double d = 0;   //c*eps + d
    
    int num_previous_solutions = pgd_solutions.solutions.size() - 1;
    assert(num_previous_solutions >= 0);
    MeshFunctionSharedPtr<double> solutionR = pgd_solutions.solutions.back();
    
    b += calc_integral_u_f(solutionR, definition.SOURCE_TERM);
    for(int i = 0; i < num_previous_solutions; i++)
    {
        a -= calc_integral_grad_u_grad_v(solutionR, pgd_solutions.solutions.at(i), definition.labels_full);

        b -= perms.EPS_AIR * calc_integral_grad_u_grad_v(solutionR, pgd_solutions.solutions.at(i), definition.labels_air);
        b -= perms.EPS_KARTIT * calc_integral_grad_u_grad_v(solutionR, pgd_solutions.solutions.at(i), definition.labels_kartit);
        b -= perms.EPS_EMPTY * calc_integral_grad_u_grad_v(solutionR, pgd_solutions.solutions.at(i), definition.labels_empty);
    }

    c += calc_integral_grad_u_grad_v(solutionR, solutionR, definition.labels_full);

    d += perms.EPS_AIR * calc_integral_grad_u_grad_v(solutionR, solutionR, definition.labels_air);
    d += perms.EPS_KARTIT * calc_integral_grad_u_grad_v(solutionR, solutionR, definition.labels_kartit);
    d += perms.EPS_EMPTY * calc_integral_grad_u_grad_v(solutionR, solutionR, definition.labels_empty);

    Function1D lastIterParam = pgd_solutions.parameters.back();
    Function1D newParam(lastIterParam.bound_lo, lastIterParam.bound_hi, lastIterParam.n_intervals);
    pgd_solutions.parameters.pop_back();
    for(int i = 0; i < newParam.n_points; i++)
    {
        newParam.values[i] = (a*newParam.points[i] + b) / (c*newParam.points[i] + d);
    }
    pgd_solutions.parameters.push_back(newParam);
}

int main(int argc, char* argv[])
{

    Function1D fn(0., 1., 2);
    fn.values[0] = 1;
    fn.values[1] = 2;
    fn.values[2] = 1;
    fn.print();
    std::cout << fn.int_F_F() << std::endl;
    std::cout << fn.value(0.9) << std::endl;
    Hermes::vector<int> profile_coarse(0,0,16,16,16,16,0,0);

    int active_electrode = 4;
    double eps_rel_material = 20;

    ProblemConfiguration configuration(profile_coarse, active_electrode);
    StandardPerms perms(eps_rel_material);
    ProblemDefinition_1 definition(configuration);

    // first solve for homogeneous BC only!
    definition.POTENTIAL = 0.0;
    definition.SOURCE_TERM = 1.0;


    // Load the mesh.
    MeshSharedPtr mesh(new Mesh);
    Hermes::Hermes2D::MeshReaderH2DXML mloader;
    mloader.load(definition.mesh_name, mesh);

    // Refine all elements, do it INIT_REF_NUM-times.
    for (unsigned int i = 0; i < INIT_REF_NUM; i++)
        mesh->refine_all_elements();

    MeshFunctionSharedPtr<double> sln = solve_problem(definition, perms, mesh);
    MeshFunctionSharedPtr<double> sln_perm = solve_permitivity(definition, perms, mesh);


    std::cout << "Integral is " << calc_integral_energy(sln, definition, perms) << std::endl;

    // Visualize the solution.
    Hermes::Hermes2D::Views::ScalarView viewS("Solution", new Hermes::Hermes2D::Views::WinGeom(0, 0, 1500, 700));
    Hermes::Hermes2D::Views::ScalarView viewP("Permitivity", new Hermes::Hermes2D::Views::WinGeom(0, 700, 1500, 700));
    //    Hermes::Hermes2D::Views::OrderView viewSp("Space", new Hermes::Hermes2D::Views::WinGeom(0, 600, 1200, 600));
    //    viewSp.show(space);
    viewP.show(sln_perm);
    viewS.show(sln);
    viewS.wait_for_close();
    return 0;
}
