#include "problem.h"
#include "definitions.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;

const int P_INIT = 4;                     // Uniform polynomial degree of mesh elements.
const int INIT_REF_NUM = 3;               // Number of initial uniform mesh refinements.

const double MIN_EPS = 1 * EPS0;
const double MAX_EPS = 10 * EPS0;


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

Function1D PGDSolutions::iteration_update_parameter()
{
    assert(solutions.size() == parameters.size());
    double a = 0;
    double b = 0;   //a*eps + b
    double c = 0;   //---------
    double d = 0;   //c*eps + d
    
    try
    {
        b += calc_integral_u_f(actual_solution, definition.SOURCE_TERM);
        for(int i = 0; i < solutions.size(); i++)
        {
            a -= calc_integral_grad_u_grad_v(actual_solution, solutions.at(i), definition.labels_full);

            b -= perms.EPS_AIR * calc_integral_grad_u_grad_v(actual_solution, solutions.at(i), definition.labels_air);
            b -= perms.EPS_KARTIT * calc_integral_grad_u_grad_v(actual_solution, solutions.at(i), definition.labels_kartit);
            b -= perms.EPS_EMPTY * calc_integral_grad_u_grad_v(actual_solution, solutions.at(i), definition.labels_empty);
        }

        c += calc_integral_grad_u_grad_v(actual_solution, actual_solution, definition.labels_full);

        d += perms.EPS_AIR * calc_integral_grad_u_grad_v(actual_solution, actual_solution, definition.labels_air);
        d += perms.EPS_KARTIT * calc_integral_grad_u_grad_v(actual_solution, actual_solution, definition.labels_kartit);
        d += perms.EPS_EMPTY * calc_integral_grad_u_grad_v(actual_solution, actual_solution, definition.labels_empty);
    }
    catch (Exceptions::Exception& e)
    {
        std::cout << e.info();
    }
    catch (std::exception& e)
    {
        std::cout << e.what();
    }

    Function1D newParam(actual_parameter.bound_lo, actual_parameter.bound_hi, actual_parameter.n_intervals);
    std::cout <<"(a,b,c,d): " << a << ", " << b << ", " << c << ", " << d << std::endl;
    for(int i = 0; i < newParam.n_points; i++)
    {
        newParam.values[i] = (a*newParam.points[i] + b) / (c*newParam.points[i] + d);
    }
    return newParam;
}

MeshFunctionSharedPtr<double> PGDSolutions::iteration_update_solution()
{
    WeakFormChangingPermInFull wf(this);

    // Initialize essential boundary conditions.
    Hermes::Hermes2D::DefaultEssentialBCConst<double> bc_essential_ground(definition.bc_labels_ground, 0);
    Hermes::Hermes2D::DefaultEssentialBCConst<double> bc_essential_potential(definition.bc_labels_potential, definition.POTENTIAL);
    Hermes::Hermes2D::EssentialBCs<double> bcs(Hermes::vector<EssentialBoundaryCondition<double> *> (&bc_essential_ground, &bc_essential_potential));

    SpaceSharedPtr<double> space(new Hermes::Hermes2D::H1Space<double>(mesh, &bcs, P_INIT));
    std::cout << "Ndofs: " << space->get_num_dofs() << std::endl;

    MeshFunctionSharedPtr<double> sln(new Solution<double>);
    Hermes::Hermes2D::LinearSolver<double> linear_solver(&wf, space);

    try
    {
        linear_solver.solve();
        double* sln_vector = linear_solver.get_sln_vector();

        Hermes::Hermes2D::Solution<double>::vector_to_solution(sln_vector, space, sln);
    }
    catch (Exceptions::Exception& e)
    {
        std::cout << e.info();
    }
    catch (std::exception& e)
    {
        std::cout << e.what();
    }

    return sln;
}

Hermes::Hermes2D::Views::ScalarView viewS("Solution", new Hermes::Hermes2D::Views::WinGeom(0, 0, 1500, 700));

void PGDSolutions::find_new_pair()
{
    Function1D func_init(MIN_EPS, MAX_EPS, 20, 1);
    actual_parameter = func_init;

    for(int i = 0; i < 10; i++)
    {
        MeshFunctionSharedPtr<double> new_solution = iteration_update_solution();
        actual_solution = new_solution;
        viewS.show(actual_solution);
        Function1D new_parameter = iteration_update_parameter();
        std::cout << "iteration " << i << ", difference " << new_parameter.diference(actual_parameter) << std::endl;
        new_parameter.print_short();
        actual_parameter = new_parameter;
    }
    parameters.push_back(actual_parameter);
    solutions.push_back(actual_solution);

}

PGDSolutions pgd_run(ProblemDefinition definition, Perms perms, MeshSharedPtr mesh)
{

    PGDSolutions pgd_solutions(definition, perms, mesh);
    for(int i = 0; i < 5; i++)
    {
        pgd_solutions.find_new_pair();

        std::cout << "NUMBER OF SOLS " << pgd_solutions.solutions.size() << std::endl;

        viewS.save_numbered_screenshot("pic/solution%03d.png", i);
        Function1D last_param = pgd_solutions.parameters.back();

        char filename[30];
        sprintf(filename, "pic/parameter%03d.dat", i);
        FILE* file;
        file = fopen(filename, "w");
        for(int i = 0; i < last_param.n_points; i++)
        {
            fprintf(file, "%g  %g\n", last_param.points[i], last_param.values[i]);
        }
        fclose(file);
    }
    return pgd_solutions;
}

const double test_x = 0.2;
const double test_y = 0.3;

void pgd_results(PGDSolutions pgd_solutions)
{
    ProblemDefinition definition = pgd_solutions.definition;
    Perms perms = pgd_solutions.perms;
    MeshSharedPtr mesh = pgd_solutions.mesh;

    FILE* file_norm = fopen("pic/normal_calculations.dat", "w");
    FILE* file_pgd = fopen("pic/pgd_calculations.dat", "w");

    for(double eps = MIN_EPS; eps <= MAX_EPS; eps += (MAX_EPS - MIN_EPS) / 10)
    {
        perms.EPS_FULL = eps;
        MeshFunctionSharedPtr<double> ref_sln = solve_problem(definition, perms, mesh);
        double val_ref = ref_sln->get_pt_value(test_x, test_y)->val[0];
        fprintf(file_norm, "%g  %g\n", eps/EPS0, val_ref);
    }
    for(double eps = MIN_EPS; eps <= MAX_EPS; eps += (MAX_EPS - MIN_EPS) / 30)
    {
        fprintf(file_pgd, "%g  ", eps/EPS0);
        for(int num_modes = 1; num_modes <= pgd_solutions.solutions.size(); num_modes++)
        {
            double val = pgd_solutions.get_pt_value(test_x, test_y, eps, num_modes);
            fprintf(file_pgd, "%g  ", val);
        }
        fprintf(file_pgd, "\n");
    }

    fclose(file_norm);
    fclose(file_pgd);
}

void simple_run(ProblemDefinition definition, Perms perms, MeshSharedPtr mesh)
{
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
    Function1D::test();
    //assert(0);

    Hermes::vector<int> profile_coarse(0,0,16,16,16,16,0,0);

    int active_electrode = 4;
    double eps_rel_material = 7;

    ProblemConfiguration configuration(profile_coarse, active_electrode);
    StandardPerms perms(eps_rel_material);
    //ProblemDefinition_Cidlo_1 definition(configuration);
    ProblemDefinition_Unit_Square definition(configuration);

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

    //simple_run(definition, perms, mesh);
    PGDSolutions pgd_solutions = pgd_run(definition, perms, mesh);
    pgd_results(pgd_solutions);

    return 0;
}
