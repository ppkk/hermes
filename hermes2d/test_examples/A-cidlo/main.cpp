#include "problem.h"
#include "definitions.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;

const double MIN_EPS = 1 * EPS0;
const double MAX_EPS = 10 * EPS0;

const int NUM_STEPS = 5;
const int MAX_STEP_ITERATIONS = 20;
const double STEP_ITERATIONS_TOLERANCE = 1e-1;

const double test_x = 0.03;
const double test_y = 0.02;


MeshFunctionSharedPtr<double> solve_problem(ProblemDefinition definition, Perms perms, MeshSharedPtr mesh, bool external_dirichlet_lift)
{
    // Initialize essential boundary conditions.
    Hermes::Hermes2D::DefaultEssentialBCConst<double> bc_essential_ground(definition.bc_labels_ground, 0);
    double potential = external_dirichlet_lift ? 0 : definition.POTENTIAL;
    Hermes::Hermes2D::DefaultEssentialBCConst<double> bc_essential_potential(definition.bc_labels_potential, potential);
    Hermes::Hermes2D::EssentialBCs<double> bcs(Hermes::vector<EssentialBoundaryCondition<double> *> (&bc_essential_ground, &bc_essential_potential));

    SpaceSharedPtr<double> space(new Hermes::Hermes2D::H1Space<double>(mesh, &bcs, definition.P_INIT));
    std::cout << "Ndofs: " << space->get_num_dofs() << std::endl;

    CustomWeakFormPoisson wf(definition, perms, external_dirichlet_lift);
    MeshFunctionSharedPtr<double> dirichlet_lift;
    if(external_dirichlet_lift)
    {
        dirichlet_lift = solve_problem(definition, AllOnePerms(), mesh, false);
        wf.set_ext(dirichlet_lift);
    }

    MeshFunctionSharedPtr<double> sln(new Solution<double>);
    Hermes::Hermes2D::LinearSolver<double> linear_solver(&wf, space);

    try
    {
        linear_solver.solve();
        double* sln_vector = linear_solver.get_sln_vector();

        Hermes::Hermes2D::Solution<double>::vector_to_solution(sln_vector, space, sln);
        if(external_dirichlet_lift)
        {

            SumFilter<double>* sum = new SumFilter<double>(Hermes::vector<MeshFunctionSharedPtr<double> >(sln, dirichlet_lift));
            return MeshFunctionSharedPtr<double>(sum);
        }
        else
        {
            return sln;
        }
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
    double a[solutions.size()], b[solutions.size()];
    double c, d, r3;
    
    try
    {
        r3 = calc_integral_u_f(actual_solution, definition.SOURCE_TERM);
        for(int i = 0; i < solutions.size(); i++)
        {
            a[i] = calc_integral_grad_u_grad_v(actual_solution, solutions.at(i), definition.labels_full);

            b[i]  = perms.EPS_AIR * calc_integral_grad_u_grad_v(actual_solution, solutions.at(i), definition.labels_air);
            b[i] += perms.EPS_KARTIT * calc_integral_grad_u_grad_v(actual_solution, solutions.at(i), definition.labels_kartit);
            b[i] += perms.EPS_EMPTY * calc_integral_grad_u_grad_v(actual_solution, solutions.at(i), definition.labels_empty);
        }

        c = calc_integral_grad_u_grad_v(actual_solution, actual_solution, definition.labels_full);

        d  = perms.EPS_AIR * calc_integral_grad_u_grad_v(actual_solution, actual_solution, definition.labels_air);
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
    for(int point_idx = 0; point_idx < newParam.n_points; point_idx++)
    {
        double epsilon = newParam.points[point_idx];

        double nominator = r3;
        for(int previous_sol = 0; previous_sol < solutions.size(); previous_sol++)
        {
            nominator -= (a[previous_sol] * epsilon +  b[previous_sol]) * parameters[previous_sol].value(epsilon);
        }
        newParam.values[point_idx] = nominator / (c * epsilon + d);
    }

    // !!!!!!!
    // delim cely vektor prvnim clenem
    // tak aby vzdy zacinal jednickou
    // zrejme by se nemelo delat az bude nehomogenni Dirichletova podminka
    //  !!!!!!!
    //newParam.normalize_first_to_one();

    return newParam;
}

MeshFunctionSharedPtr<double> PGDSolutions::iteration_update_solution()
{
    WeakFormChangingPermInFull wf(this);

    // Initialize essential boundary conditions.
    Hermes::Hermes2D::DefaultEssentialBCConst<double> bc_essential_ground(definition.bc_labels_ground, 0);
    Hermes::Hermes2D::DefaultEssentialBCConst<double> bc_essential_potential(definition.bc_labels_potential, definition.POTENTIAL);
    Hermes::Hermes2D::EssentialBCs<double> bcs(Hermes::vector<EssentialBoundaryCondition<double> *> (&bc_essential_ground, &bc_essential_potential));

    SpaceSharedPtr<double> space(new Hermes::Hermes2D::H1Space<double>(mesh, &bcs, definition.P_INIT));
    std::cout << "Ndofs: " << space->get_num_dofs() << std::endl;

    MeshFunctionSharedPtr<double> sln(new Solution<double>);
    Hermes::Hermes2D::LinearSolver<double> linear_solver(&wf, space);

    try
    {
        linear_solver.solve();
        double* sln_vector = linear_solver.get_sln_vector();
        double max = 0;
        for(int i = 0; i < space->get_num_dofs(); i++)
            max = std::max(max, abs(sln_vector[i]));
        printf("maximal sol comp %g\n", max);

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

FILE* convergence_file;

Hermes::Hermes2D::Views::ScalarView viewS("Solution", new Hermes::Hermes2D::Views::WinGeom(0, 0, 1500, 700));

void PGDSolutions::find_new_pair()
{
    Function1D func_init(MIN_EPS, MAX_EPS, 100, 1);
    actual_parameter = func_init;

    double difference = 1e10;
    int iteration = 1;
    while((difference > STEP_ITERATIONS_TOLERANCE) && (iteration < MAX_STEP_ITERATIONS))
    {
        std::cout << "iteration " << iteration << std::endl;
        MeshFunctionSharedPtr<double> new_solution = iteration_update_solution();
        actual_solution = new_solution;
        viewS.show(actual_solution);
        Function1D new_parameter = iteration_update_parameter();
        difference = new_parameter.diference(actual_parameter) / new_parameter.norm();
        fprintf(convergence_file, "%g ", difference);
        //new_parameter.print_short();
        actual_parameter = new_parameter;
        iteration++;
    }
    fprintf(convergence_file, "\n");
    parameters.push_back(actual_parameter);
    solutions.push_back(actual_solution);

}

PGDSolutions pgd_run(ProblemDefinition definition, Perms perms, MeshSharedPtr mesh)
{
    convergence_file = fopen("pic/convergence.dat", "w");

    PGDSolutions pgd_solutions(definition, perms, mesh);
    for(int i = 0; i < NUM_STEPS; i++)
    {
        pgd_solutions.find_new_pair();

        std::cout << "NUMBER OF SOLS " << pgd_solutions.solutions.size() << std::endl;

        viewS.save_numbered_screenshot("pic/sol_mode%03d.png", i);
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

    fclose(convergence_file);
    return pgd_solutions;
}

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
        MeshFunctionSharedPtr<double> ref_sln = solve_problem(definition, perms, mesh, false);
        double val_ref = ref_sln->get_pt_value(test_x, test_y)->val[0];
        fprintf(file_norm, "%g  %g\n", eps/EPS0, val_ref);

        viewS.show(ref_sln);
        char sol_name[20];
        sprintf(sol_name, "pic/solution_perm_%1.1lf.png", int(10*eps/EPS0)/10.);
        viewS.save_screenshot(sol_name);

        MeshFunctionSharedPtr<double>  pgd_sln = pgd_solutions.get_filter(eps);
        viewS.show(pgd_sln);
        sprintf(sol_name, "pic/pgd_solution_perm_%1.1lf.png", int(10*eps/EPS0)/10.);
        viewS.save_screenshot(sol_name);
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

void simple_run(ProblemDefinition definition, Perms perms, MeshSharedPtr mesh, bool external_dirichlet_lift)
{
    MeshFunctionSharedPtr<double> sln = solve_problem(definition, perms, mesh, external_dirichlet_lift);
    MeshFunctionSharedPtr<double> sln_perm = solve_permitivity(definition, perms, mesh);

    //std::cout << "Integral is " << calc_integral_energy(sln, definition, perms) << std::endl;

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
    Function1D::test();

    Hermes::vector<int> profile_coarse(0,0,16,16,16,16,0,0);

    int active_electrode = 4;
    double eps_rel_material = 7;

    ProblemConfiguration configuration(profile_coarse, active_electrode);
    StandardPerms perms(eps_rel_material);
    ProblemDefinitionCidlo1 definition(configuration);
    //ProblemDefinitionUnitSquare definition(configuration);
    //ProblemDefinitionUnitSquareDivided definition(configuration);

    // first solve for homogeneous BC only!
    //definition.POTENTIAL = 0.0;
    //definition.SOURCE_TERM = 1.0;


    // Load the mesh.
    MeshSharedPtr mesh(new Mesh);
    Hermes::Hermes2D::MeshReaderH2DXML mloader;
    mloader.load(definition.mesh_name, mesh);

    // Refine all elements, do it INIT_REF_NUM-times.
    for (unsigned int i = 0; i < definition.INIT_REF_NUM; i++)
        mesh->refine_all_elements();

    simple_run(definition, perms, mesh, true);
    //PGDSolutions pgd_solutions = pgd_run(definition, perms, mesh);
    //pgd_results(pgd_solutions);

    return 0;
}
