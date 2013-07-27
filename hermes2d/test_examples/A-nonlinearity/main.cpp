#include <ctime>

#include "definitions.h"

//  This example solves the same nonlinear problem as the previous
//  one but now using the Newton's method.
//
//  PDE: Stationary heat transfer equation with nonlinear thermal
//       conductivity, - div[lambda(u) grad u] + src(x, y) = 0.
//
//  Nonlinearity: lambda(u) = 1 + Hermes::pow(u, alpha).
//
//  Domain: square (-10, 10)^2.
//
//  BC: Nonconstant Dirichlet.
//
//  The following parameters can be changed:

// Initial polynomial degree.
const int P_INIT = 2;                             
// Stopping criterion for the Newton's method.
const double NEWTON_TOL = 1e-8;                   
// Maximum allowed number of Newton iterations.
const int NEWTON_MAX_ITER = 100;                  
// Number of initial uniform mesh refinements.
const int INIT_GLOB_REF_NUM = 3;                  
// Number of initial refinements towards boundary.
const int INIT_BDY_REF_NUM = 4;                   

// Problem parameters.
double heat_src = 1.0;
double alpha = 18.0;

template <typename Scalar>
class MyNewton : public Hermes::Hermes2D::NewtonSolver<Scalar>
{
public:
    MyNewton(DiscreteProblem<Scalar>* dp) : Hermes::Hermes2D::NewtonSolver<Scalar>(dp) {}

    virtual bool on_finish();

    int iterations;
    int jacobianRecalcs;
};

template <typename Scalar>
bool MyNewton<Scalar>::on_finish()
{
    const Hermes::vector<bool>& jacobian_recalculated_log = this->get_parameter_value(this->iterations_with_recalculated_jacobian());

    iterations = jacobian_recalculated_log.size();
    //newton.get_parameter_value(newton.iteration());
    int iteration2 = this->get_current_iteration_number();
    assert(iterations + 1 == iteration2);
    //std::cout << "\n\n AAAAAA " << iterations << ", " << iteration2 << "\n\n\n";

    jacobianRecalcs = 0;
    for(int ijac = 0; ijac < jacobian_recalculated_log.size(); ijac++)
        if(jacobian_recalculated_log[ijac])
            jacobianRecalcs++;
    //jacobianRecalcs = newton.get_parameter_value(newton.iteration_with_recalculated_jacobian());
}


int main(int argc, char* argv[])
{
#pragma region initialization

    // Load the mesh.
    MeshSharedPtr mesh(new Mesh);
    MeshReaderH2D mloader;
    mloader.load("square.mesh", mesh);

    // Perform initial mesh refinements.
    for(int i = 0; i < INIT_GLOB_REF_NUM; i++) mesh->refine_all_elements();
    mesh->refine_towards_boundary("Bdy", INIT_BDY_REF_NUM);

    // Initialize boundary conditions.
    CustomEssentialBCNonConst bc_essential("Bdy");
    EssentialBCs<double> bcs(&bc_essential);

    // Create an H1 space with default shapeset.
    SpaceSharedPtr<double> space(new H1Space<double>(mesh, &bcs, P_INIT));
    int ndof = space->get_num_dofs();
    Hermes::Mixins::Loggable::Static::info("ndof: %d", ndof);

    // Initialize the weak formulation
    CustomNonlinearity lambda(alpha);
    Hermes2DFunction<double> src(-heat_src);
    DefaultWeakFormPoisson<double> wf(HERMES_ANY, &lambda, &src);


    FILE *file;
    file = fopen("results.txt", "w");

//    double suff_imp_jac_step = 0.05;
//    double suff_imp_jac_num = 20;
//    double suff_imp_jac_min = 0.05;

    double suff_imp_jac = 0.001;

    double suff_imp_step = 0.05;
    double suff_imp_num = 18;
    double suff_imp_min = 0.9;

    double suff_imp = 1;

    double initial_damp_step = 0.05;
    double initial_damp_num = 16;
    double initial_damp_min = 0.25;

    double initial_damp = 1;

    double max_st_reuse = 0;


//    suff_imp_jac = suff_imp_jac_min;
//    for(int suff_imp_jac_idx = 0; suff_imp_jac_idx < suff_imp_jac_num; suff_imp_jac_idx++)
    suff_imp = suff_imp_min;
    for(int suff_imp_idx = 0; suff_imp_idx < suff_imp_num; suff_imp_idx++)
    {
        initial_damp = initial_damp_min;
        for(int initial_damp_idx = 0; initial_damp_idx < initial_damp_num; initial_damp_idx++)
        {
            clock_t begin = clock();


            // Initialize the FE problem.
            DiscreteProblem<double> dp(&wf, space);

            // Project the initial condition on the FE space to obtain initial
            // coefficient vector for the Newton's method.
            // NOTE: If you want to start from the zero vector, just define
            // coeff_vec to be a vector of ndof zeros (no projection is needed).
            Hermes::Mixins::Loggable::Static::info("Projecting to obtain initial vector for the Newton's method.");
            double* coeff_vec = new double[ndof];
            MeshFunctionSharedPtr<double> init_sln(new CustomInitialCondition(mesh));
            OGProjection<double> ogProjection; ogProjection.project_global(space, init_sln, coeff_vec);

#pragma endregion

            // Initialize Newton solver.
            MyNewton<double> newton(&dp);
            newton.set_convergence_measurement(SolutionDistanceFromPreviousRelative | ResidualNormAbsolute);
            newton.set_tolerance(1e-2, SolutionDistanceFromPreviousRelative, true);
            newton.set_tolerance(1e-3, ResidualNormAbsolute, true);
            newton.set_min_allowed_damping_coeff(1e-12);
            newton.set_max_allowed_residual_norm(1e159);
            newton.set_max_allowed_iterations(1000);

            // muj pokus
            newton.set_sufficient_improvement_factor_jacobian(suff_imp_jac);
            newton.set_max_steps_with_reused_jacobian(max_st_reuse);
            newton.set_initial_auto_damping_coeff(initial_damp);
            newton.set_sufficient_improvement_factor(suff_imp);


            // 1st - OK
            //  newton.set_sufficient_improvement_factor_jacobian(0.05);
            //  newton.set_max_steps_with_reused_jacobian(99);
            //  newton.set_initial_auto_damping_coeff(.5);
            //  newton.set_sufficient_improvement_factor(1.25);

            // 2nd - OK
            // newton.set_sufficient_improvement_factor_jacobian(0.9);
            // newton.set_max_steps_with_reused_jacobian(10);
            // newton.set_manual_damping_coeff(0.1);

            bool converged;

            converged = true;
            // Perform Newton's iteration.
            try
            {
                newton.solve(coeff_vec);
            }
            catch(NewtonSolver<double>::NewtonException& e)
            {
                converged = false;
                switch(e.get_exception_state())
                {
                case NewtonSolver<double>::AboveMaxIterations:
                    std::cout << std::endl << "\t\t\tAboveMaxIterations" << std::endl;
                    break;
                case NewtonSolver<double>::BelowMinDampingCoeff:
                    std::cout << std::endl << "\t\t\tBelowMinDampingCoeff" << std::endl;
                    break;
                case NewtonSolver<double>::AboveMaxAllowedResidualNorm:
                    std::cout << std::endl << "\t\t\tAboveMaxAllowedResidualNorm" << std::endl;
                    break;
                }
            }

#pragma region finalization

            catch(std::exception& e)
            {
                converged = false;
                std::cout << e.what();
            }

            int iteration = 1000;
            int jacobianRecalcs = 1000;

            if(converged)
            {
                iteration = newton.iterations;
                jacobianRecalcs = newton.jacobianRecalcs;
                if(iteration == 1)
                {
                    iteration = 1000;
                    jacobianRecalcs = 1000;
                }
            }

            clock_t end = clock();
            double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;


            fprintf(file, "%lf %lf %lf %lf %d %d %lf\n",  suff_imp_jac, max_st_reuse, initial_damp, suff_imp, iteration, jacobianRecalcs, elapsed_secs);

           // Clean up.
           delete [] coeff_vec;


           initial_damp += initial_damp_step;
           if(initial_damp > 1.0)
               initial_damp = 1.0;
        }
        suff_imp += suff_imp_step;
    }

    fclose(file);

//    // Translate the resulting coefficient vector into a Solution.
//    MeshFunctionSharedPtr<double> sln(new Solution<double>);
//    Solution<double>::vector_to_solution(newton.get_sln_vector(), space, sln);

//    // Visualise the solution and mesh.
//    ScalarView s_view("Solution", new WinGeom(0, 0, 440, 350));
//    s_view.show_mesh(false);
//    s_view.show(sln);
//    OrderView o_view("Mesh", new WinGeom(450, 0, 400, 350));
//    o_view.show(space);

//    // Wait for all views to be closed.
//    View::wait();
//    return 0;

#pragma endregion
}
