//#include "solvers/interfaces/umfpack_solver.h"

#include "definitions.h"
#include "cg.h"
#include "stdio.h"
#include "gmres.h"
#include "iml_wrappers.h"

using namespace RefinementSelectors;

//  This is another example that allows you to compare h- and hp-adaptivity from the point of view
//  of both CPU time requirements and discrete problem size, look at the quality of the a-posteriori
//  error estimator used by Hermes (exact error is provided), etc. You can also change
//  the parameter MESH_REGULARITY to see the influence of hanging nodes on the adaptive process.
//  The problem is made harder for adaptive algorithms by increasing the parameter SLOPE.
//
//  PDE: -Laplace u + f = 0
//
//  Known exact solution.
//
//  Domain: unit square (0, 1) x (0, 1), see the file square.mesh->
//
//  BC:  Dirichlet, given by exact solution.
//
//  The following parameters can be changed:

// Initial polynomial degree of mesh elements.
const int P_INIT = 1;
// Number of initial uniform mesh refinements.
const int INIT_REF_NUM = 2;

const int REFERENCE_ORDER_INCREASE = 0;
// This is a quantitative parameter of the adapt(...) function and
// it has different meanings for various adaptive strategies.
const double THRESHOLD = 0.9;
// Error calculation & adaptivity.
DefaultErrorCalculator<double, HERMES_H1_NORM> errorCalculator(RelativeErrorToGlobalNorm, 1);
// Stopping criterion for an adaptivity step.
AdaptStoppingCriterionSingleElement<double> stoppingCriterion(THRESHOLD);
// Adaptivity processor class.
Adapt<double> adaptivity(&errorCalculator, &stoppingCriterion);
// Predefined list of element refinement candidates.
//const CandList CAND_LIST = H2D_HP_ANISO;
const CandList CAND_LIST = H2D_H_ISO;
// Stopping criterion for adaptivity.
const double ERR_STOP = 1e-1;

// interpolate or project?
bool interpolate = true;

// Problem parameters.
// Slope of the layer.
double slope = 60;                                

int file_index = 0;
double minr = -2, maxr = 2;
Views::LinearizerCriterionFixed lincrit(4);
SimpleVector<double> last_directly_solved_rhs;

class ProjectionPrecond
{
public:
    ProjectionPrecond()
    {
        show_progress = false;
    }
    IMLVector solve(const IMLVector& x) const
    {
        SimpleVector<double> vector_in;
        x.copy(&vector_in);

        Views::OrderView oview_coarse("Coarse orders", new Views::WinGeom(600, 0, 420, 350));
        Views::OrderView oview_fine("Fine orders", new Views::WinGeom(1100, 0, 420, 350));
        Views::ScalarView sview_fine("Fine", new Views::WinGeom(0, 400, 440, 350));
        sview_fine.show_mesh(true);
        sview_fine.fix_scale_width(50);
        Views::ScalarView sview_coarse_projection("Coarse projection", new Views::WinGeom(0, 800, 440, 350));
        sview_coarse_projection.show_mesh(true);
        sview_coarse_projection.fix_scale_width(50);
        Views::ScalarView sview_coarse_solution("Coarse solution", new Views::WinGeom(500, 800, 440, 350));
        sview_coarse_solution.show_mesh(true);
        sview_coarse_solution.fix_scale_width(50);
        Views::ScalarView sview_back("Back", new Views::WinGeom(500, 400, 440, 350));
        sview_back.show_mesh(true);
        sview_back.fix_scale_width(50);
        Views::ScalarView sview_exact_precond("Exact precond", new Views::WinGeom(950, 400, 440, 350));
        sview_exact_precond.show_mesh(true);
        sview_exact_precond.fix_scale_width(50);
        Views::ScalarView sview_other("Other", new Views::WinGeom(950, 800, 440, 350));
        sview_other.show_mesh(true);
        sview_other.fix_scale_width(50);


        sview_fine.set_min_max_range(minr,maxr);
        sview_coarse_projection.set_min_max_range(minr,maxr);
        sview_coarse_solution.set_min_max_range(minr,maxr);
        sview_back.set_min_max_range(minr,maxr);
        sview_exact_precond.set_min_max_range(minr,maxr);
        sview_other.set_min_max_range(minr,maxr);

        sview_fine.get_linearizer()->set_criterion(lincrit);
        sview_coarse_projection.get_linearizer()->set_criterion(lincrit);
        sview_coarse_solution.get_linearizer()->set_criterion(lincrit);
        sview_back.get_linearizer()->set_criterion(lincrit);
        sview_exact_precond.get_linearizer()->set_criterion(lincrit);
        sview_other.get_linearizer()->set_criterion(lincrit);

        // fine direct solver

        char file_name[100];
//        sprintf(file_name, "vector_in_%d", file_index);
//        vector_in.export_to_file(file_name, "in", EXPORT_FORMAT_MATRIX_MARKET, "%g");
        MeshFunctionSharedPtr<double> tmp_solution(new Solution<double>()), tmp_solution2(new Solution<double>());
        Solution<double>::vector_to_solution(vector_in.v, fine_space, tmp_solution);
        if(show_progress) sview_fine.show(tmp_solution);

        SimpleVector<double> projection;
        projection.alloc(coarse_space->get_num_dofs());
        OGProjection<double>::project_global(coarse_space, tmp_solution, &projection);

        if(show_progress) Solution<double>::vector_to_solution(projection.v, coarse_space, tmp_solution);
        if(show_progress) sview_coarse_projection.show(tmp_solution);

        rhs->alloc(projection.get_size());
        rhs->set_vector(&projection);
        // zkouska: misto rezidualu pouzit presnou hodnotu prave strany pro posledni primo reseny problem
        if(show_progress) Solution<double>::vector_to_solution(last_directly_solved_rhs.v, coarse_space, tmp_solution);
        if(show_progress) sview_other.show(tmp_solution);
//        rhs->alloc(last_directly_solved_rhs.get_size());
//        rhs->set_vector(&last_directly_solved_rhs);

        direct_solver->set_reuse_scheme(Hermes::Solvers::HERMES_CREATE_STRUCTURE_FROM_SCRATCH);
        direct_solver->solve();
        Solution<double>::vector_to_solution(direct_solver->get_sln_vector(), coarse_space, tmp_solution2);

        if(show_progress) sview_coarse_solution.show(tmp_solution2);

        SimpleVector<double> result;
        result.alloc(fine_space->get_num_dofs());
        OGProjection<double>::project_global(fine_space, tmp_solution2, &result);


        // comparison with exact preconditioner
        rhs_fine->set_vector(&vector_in);
        direct_solver_fine->set_reuse_scheme(Hermes::Solvers::HERMES_CREATE_STRUCTURE_FROM_SCRATCH);
        direct_solver_fine->solve();

        if(show_progress) Solution<double>::vector_to_solution(direct_solver_fine->get_sln_vector(), fine_space, tmp_solution);
        if(show_progress) sview_exact_precond.show(tmp_solution);
        //result.set_vector(direct_solver_fine->get_sln_vector());


//        // or finer solution moved to coarse mesh and back
//        projection.alloc(coarse_space->get_num_dofs());
//        Solution<double>::vector_to_solution(result.v, fine_space, tmp_solution);
//        OGProjection<double>::project_global(coarse_space, tmp_solution, &projection);
//        result.alloc(fine_space->get_num_dofs());
//        Solution<double>::vector_to_solution(projection.v, coarse_space, tmp_solution2);
//        OGProjection<double>::project_global(fine_space, tmp_solution2, &result);

//        // or try to exactly solve this perturbed projection
//        rhs_fine->set_vector(vector_in); //result);
//        direct_solver_fine->set_reuse_scheme(Hermes::Solvers::HERMES_CREATE_STRUCTURE_FROM_SCRATCH);
//        direct_solver_fine->solve();
//        result->set_vector(direct_solver_fine->get_sln_vector());


        // or not to do anything at all:
//        result.set_vector(&vector_in);

        if(show_progress) Solution<double>::vector_to_solution(result.v, fine_space, tmp_solution);
        if(show_progress) sview_back.show(tmp_solution);


//        sprintf(file_name, "vector_out_%d", file_index++);
//        result.export_to_file(file_name, "out", EXPORT_FORMAT_MATRIX_MARKET, "%g");

        if(show_progress) oview_coarse.show(coarse_space);
        if(show_progress) oview_fine.show(fine_space);
        if(show_progress) getchar();

        IMLVector imlResult;
        imlResult.set(&result);
        return imlResult;
    }

    Solvers::UMFPackLinearMatrixSolver<double>* direct_solver;
    SimpleVector<double>* rhs;
    SpaceSharedPtr<double> coarse_space;

    Solvers::UMFPackLinearMatrixSolver<double>* direct_solver_fine;
    SimpleVector<double>* rhs_fine;
    SpaceSharedPtr<double> fine_space;

    bool show_progress;

};

void project_or_interpolate(SpaceSharedPtr<double> orig_space, MeshFunctionSharedPtr<double> orig_solution, SpaceSharedPtr<double> target_space, double* target_sln_vec)
{
    if(interpolate)
    {
        assert(0);
        //project_or_interpolate(orig_space, orig_solution, target_space, target_sln_vec);
    }
    else
    {
        OGProjection<double>::project_global(target_space, orig_solution, target_sln_vec);
    }
}

void project_or_interpolate(SpaceSharedPtr<double> orig_space, MeshFunctionSharedPtr<double> orig_solution, SpaceSharedPtr<double> target_space, SimpleVector<double>* target_sln_vec)
{
    if(interpolate)
    {
        assert(0);
        //orig_solution.get_solution()->get_v
        //project_or_interpolate(orig_space, orig_solution, target_space, target_sln_vec);
    }
    else
    {
        OGProjection<double>::project_global(target_space, orig_solution, target_sln_vec);
    }
}


void project_or_interpolate(SpaceSharedPtr<double> orig_space, double* orig_sln_vec, SpaceSharedPtr<double> target_space, SimpleVector<double>* target_sln_vec)
{
    if(interpolate)
    {
        VertexBasedInterpolation<double>::interpolate(orig_space, orig_sln_vec, target_space, target_sln_vec->v);
    }
    else
    {
        MeshFunctionSharedPtr<double> tmp_solution(new Solution<double>);
        Solution<double>::vector_to_solution(orig_sln_vec, orig_space, tmp_solution);
        project_or_interpolate(orig_space, tmp_solution, target_space, target_sln_vec);
    }
}

void project_or_interpolate(SpaceSharedPtr<double> orig_space, double* orig_sln_vec, SpaceSharedPtr<double> target_space, double* target_sln_vec)
{
    if(interpolate)
    {
        VertexBasedInterpolation<double>::interpolate(orig_space, orig_sln_vec, target_space, target_sln_vec);
    }
    else
    {
        MeshFunctionSharedPtr<double> tmp_solution(new Solution<double>);
        Solution<double>::vector_to_solution(orig_sln_vec, orig_space, tmp_solution);
        project_or_interpolate(orig_space, tmp_solution, target_space, target_sln_vec);
    }
}

void test_projections(SpaceSharedPtr<double> coarse_space, SpaceSharedPtr<double> fine_space)
{
    Views::ScalarView sview_original("Original", new Views::WinGeom(0, 800, 440, 350));
    sview_original.show_mesh(true);
    sview_original.fix_scale_width(50);
    Views::ScalarView sview_projection("Projection", new Views::WinGeom(500, 800, 440, 350));
    sview_projection.show_mesh(true);
    sview_projection.fix_scale_width(50);
    Views::ScalarView sview_projection_by_matrix("Projection by matrix", new Views::WinGeom(950, 800, 440, 350));
    sview_projection_by_matrix.show_mesh(true);
    sview_projection_by_matrix.fix_scale_width(50);
    sview_original.set_min_max_range(minr,maxr);
    sview_projection.set_min_max_range(minr,maxr);
    sview_projection_by_matrix.set_min_max_range(minr,maxr);
    sview_original.get_linearizer()->set_criterion(lincrit);
    sview_projection.get_linearizer()->set_criterion(lincrit);
    sview_projection_by_matrix.get_linearizer()->set_criterion(lincrit);

    int n_coarse = coarse_space->get_num_dofs();
    int n_fine = fine_space->get_num_dofs();

    IMLMatrix coarsen(n_coarse, n_fine);
    SimpleVector<double> projection;
    MeshFunctionSharedPtr<double> tmp_solution(new Solution<double>);
    projection.alloc(n_coarse);
    double* vec = new double[n_fine];
    for(int i = 0; i < n_fine; i++)
    {
        memset(vec, 0, n_fine*sizeof(double));
        vec[i] = 1.;

        projection.zero();
        project_or_interpolate(fine_space, vec, coarse_space, &projection);
        for(int j = 0; j < n_coarse; j++)
            coarsen(j, i) = projection.v[j];

        Solution<double>::vector_to_solution(vec, fine_space, tmp_solution);
        sview_original.show(tmp_solution);
        Solution<double>::vector_to_solution(projection.v, coarse_space, tmp_solution);
        sview_projection.show(tmp_solution);
        getchar();
    }

    IMLMatrix refine(n_fine, n_coarse);
    projection.alloc(n_fine);
    double* vec2 = new double[n_coarse];
    for(int i = 0; i < n_coarse; i++)
    {
        memset(vec, 0, n_coarse*sizeof(double));
        vec[i] = 1.;

        projection.zero();
        project_or_interpolate(coarse_space, vec, fine_space, &projection);
        for(int j = 0; j < n_fine; j++)
            refine(j, i) = projection.v[j];

        Solution<double>::vector_to_solution(vec, coarse_space, tmp_solution);
        sview_original.show(tmp_solution);
        Solution<double>::vector_to_solution(projection.v, fine_space, tmp_solution);
        sview_projection.show(tmp_solution);
       getchar();
    }

    projection.alloc(n_coarse);
    memset(vec, 0, n_fine*sizeof(double));
    for(int i = 0; i < n_fine; i++)
        vec[i] = i%4; //1
    IMLVector iml_vec;
    iml_vec.set(vec, n_fine);
    Solution<double>::vector_to_solution(vec, fine_space, tmp_solution);
    sview_original.show(tmp_solution);
    OGProjection<double>::project_global(coarse_space, tmp_solution, &projection);
    Solution<double>::vector_to_solution(projection.v, coarse_space, tmp_solution);
    sview_projection.show(tmp_solution);
    IMLVector iml_projection;
    iml_projection.set(&projection);
    IMLVector iml_vec_2 = coarsen * iml_vec;
    std::cout << " je projekce linearni? : " << norm(iml_vec_2 - iml_projection) << std::endl;
    Solution<double>::vector_to_solution(iml_vec_2.m_data, coarse_space, tmp_solution);
    sview_projection_by_matrix.show(tmp_solution);
    getchar();

    projection.alloc(n_fine);
    memset(vec, 0, n_coarse*sizeof(double));
    for(int i = 0; i < n_coarse; i++)
        vec[i] = i%4; //1
    iml_vec.set(vec, n_coarse);
    Solution<double>::vector_to_solution(vec, coarse_space, tmp_solution);
    sview_original.show(tmp_solution);
    OGProjection<double>::project_global(fine_space, tmp_solution, &projection);
    Solution<double>::vector_to_solution(projection.v, fine_space, tmp_solution);
    sview_projection.show(tmp_solution);
    iml_projection.set(&projection);
    iml_vec_2.alloc(n_fine);
    iml_vec_2 = refine * iml_vec;
    std::cout << " je projekce linearni? : " << norm(iml_vec_2 - iml_projection) << std::endl;
    Solution<double>::vector_to_solution(iml_vec_2.m_data, fine_space, tmp_solution);
    sview_projection_by_matrix.show(tmp_solution);
    getchar();

    coarsen.print_sparse("/home/pkus/matlab/multigrid/restrict.m", "R");
    refine.print_sparse("/home/pkus/matlab/multigrid/interpolate.m", "I");

    getchar();
}


void test_projection_precond(const ProjectionPrecond& precond, const IMLPrecondTimesMatrixOperator<IMLOperator, ProjectionPrecond> & precondOperator)
{
    int size = precond.fine_space->get_num_dofs();
    IMLVector vec(size), result1(size), result2(size);
    for(int i = 0; i < size; i++)
    {
        vec(i) = i;
    }

    result1 = precond.solve(vec);
    result2 = precond.solve(vec);

    std::cout << "test precond: norm = " << norm(result1-result2) << std::endl;

    result1 = precondOperator * vec;
    result2 = precondOperator * vec;

    std::cout << "test precond operator: norm = " << norm(result1-result2) << std::endl;
}

int main(int argc, char* argv[])
{
  //test_iml_wrappers();
  adaptivity.set_regularization_level(1);

  // Load the mesh.
  MeshSharedPtr mesh(new Mesh);
  MeshReaderH2D mloader;
  // Quadrilaterals.
  mloader.load("square_quad.mesh", mesh);
  // Triangles.
  // mloader.load("square_tri.mesh", mesh);   

  // Perform initial mesh refinements.
  for (int i = 0; i < INIT_REF_NUM; i++) mesh->refine_all_elements();
  
  // Define exact solution.
  MeshFunctionSharedPtr<double> exact_sln(new CustomExactSolution(mesh, slope));

  // Define custom function f.
  CustomFunction f(slope);

  // Initialize the weak formulation.
  //Hermes::Hermes1DFunction<double> lambda(1.0);
  DefaultWeakFormPoissonLinear<double> wf(HERMES_ANY, /*&lambda,*/ &f);
  
  // Initialize boundary conditions
  DefaultEssentialBCNonConst<double> bc_essential("Bdy", exact_sln);
  EssentialBCs<double> bcs(&bc_essential);

  // Create an H1 space with default shapeset.
  SpaceSharedPtr<double> space(new H1Space<double>(mesh, &bcs, P_INIT));
  SpaceSharedPtr<double> last_directly_solved_space(new H1Space<double>());
  MeshSharedPtr last_directly_solved_mesh(new Mesh);
  
  // Initialize approximate solution.
  MeshFunctionSharedPtr<double> sln(new Solution<double>());

  // Initialize refinement selector.
  H1ProjBasedSelector<double> selector(CAND_LIST);

  // Initialize views.
  Views::ScalarView sview_coarse("Solution coarse", new Views::WinGeom(0, 0, 440, 350));
  sview_coarse.show_mesh(false);
  sview_coarse.fix_scale_width(50);
  sview_coarse.set_min_max_range(minr,maxr);
  Views::ScalarView sview_fine("Solution fine", new Views::WinGeom(450, 0, 440, 350));
  sview_fine.show_mesh(false);
  sview_fine.fix_scale_width(50);
  sview_fine.set_min_max_range(minr,maxr);
  Views::OrderView oview("Polynomial orders", new Views::WinGeom(900, 0, 420, 350));

  // DOF and CPU convergence graphs.
  SimpleGraph graph_dof_est, graph_cpu_est, graph_dof_exact, graph_cpu_exact;

  // Time measurement.
  Hermes::Mixins::TimeMeasurable cpu_time;
  cpu_time.tick();

  CSCMatrix<double> m, m_fine;
  SimpleVector<double> rhs, rhs_fine, x;
  Solvers::UMFPackLinearMatrixSolver<double> umfpack_solver(&m, &rhs);
  Solvers::UMFPackLinearMatrixSolver<double> umfpack_solver_fine(&m_fine, &rhs_fine);

  // Adaptivity loop:
  int as = 1; bool done = false;
  do
  {
      if(as == 8) break;
    cpu_time.tick();

    // Construct globally refined reference mesh and setup reference space.
    Mesh::ReferenceMeshCreator refMeshCreator(mesh);
    MeshSharedPtr ref_mesh = refMeshCreator.create_ref_mesh();

    Space<double>::ReferenceSpaceCreator refSpaceCreator(space, ref_mesh, REFERENCE_ORDER_INCREASE);
    SpaceSharedPtr<double> ref_space = refSpaceCreator.create_ref_space();
    int ndof_ref = ref_space->get_num_dofs();

    Hermes::Mixins::Loggable::Static::info("---- Adaptivity step %d (%d DOF):", as, ndof_ref);
    cpu_time.tick();
    
    Hermes::Mixins::Loggable::Static::info("Solving on reference mesh.");

    // Assemble the discrete problem.
    DiscreteProblem<double> dp(&wf, ref_space);
    dp.set_linear();
    DiscreteProblem<double> dp_fine(&wf, ref_space);
    dp_fine.set_linear();
    Space<double>::assign_dofs(dp.get_spaces());

    //umfpack_solver.
    //double coeff_vec[ndof_ref];
    dp_fine.assemble(&m_fine, &rhs_fine);
    MeshFunctionSharedPtr<double> ref_sln(new Solution<double>());

    int gmres_m = 50;
    double tol = 1e-10;
    int max_iter = 50;
    IMLOperator iml_operator(&m_fine);
    IMLMatrix iml_matrix(gmres_m+1, ndof_ref);
    IMLVector iml_x, iml_b;
    iml_x.alloc(ndof_ref);
    iml_x = 0.0;
    iml_b.alloc(ndof_ref);
    iml_b.set(&rhs_fine);

    IMLEmptyPreconditioner empty_preconditioner;
    IMLDiagPreconditioner diag_preconditioner(&m_fine);
    ProjectionPrecond projection_preconditioner;
    projection_preconditioner.coarse_space = last_directly_solved_space;
    projection_preconditioner.fine_space = ref_space;
    projection_preconditioner.direct_solver = &umfpack_solver;
    projection_preconditioner.rhs = &rhs;
    projection_preconditioner.direct_solver_fine = &umfpack_solver_fine;
    projection_preconditioner.rhs_fine = &rhs_fine;

    IMLPrecondTimesMatrixOperator<IMLOperator, ProjectionPrecond> projection_precond_operator(&iml_operator, &projection_preconditioner);
    IMLPrecondTimesMatrixOperator<IMLOperator, IMLDiagPreconditioner> diag_precond_operator(&iml_operator, &diag_preconditioner);

    bool solve_directly = (as % 4 == 1);

    if(solve_directly)
    {
        Hermes::Mixins::Loggable::Static::info("+++++++ adaptivity step %d -> using      DIRECT     solver", as);
        dp.assemble(&m, &rhs);
        umfpack_solver.set_reuse_scheme(Hermes::Solvers::HERMES_CREATE_STRUCTURE_FROM_SCRATCH);
        umfpack_solver.solve();
        Solution<double>::vector_to_solution(umfpack_solver.get_sln_vector(), ref_space, ref_sln);
        last_directly_solved_mesh->copy(ref_space->get_mesh());
        last_directly_solved_rhs.alloc(rhs.get_size());
        last_directly_solved_rhs.set_vector(&rhs);

        IMLOperator oper(&m);
        IMLMatrix matrix(rhs.get_size(), rhs.get_size());
        IMLOperatorToMatrix(oper, matrix);
        //std::cout << "Directly solved matrix " << matrix;
        matrix.print_sparse("/home/pkus/matlab/multigrid/matrixAsmall.m", "Asmall");
        Space<double>::ReferenceSpaceCreator spaceCreator(ref_space, last_directly_solved_mesh, 0);
        last_directly_solved_space = spaceCreator.create_ref_space();
        //last_directly_solved_space->copy(ref_space, last_directly_solved_mesh);
    }
    else
    {
        int converged = 0;
        test_projections(last_directly_solved_space, ref_space);
        test_projection_precond(projection_preconditioner, projection_precond_operator);

        Hermes::Mixins::Loggable::Static::info("+++++++ adaptivity step %d -> using     ITERATIVE    solver", as);
        x.alloc(rhs_fine.get_size());
        x.zero();

        IMLMatrix dense_matrix(ndof_ref, ndof_ref);
        IMLMatrix dense_precond_matrix(ndof_ref, ndof_ref);
        IMLOperatorToMatrix(iml_operator, dense_matrix);
        dense_matrix.print_sparse("/home/pkus/matlab/multigrid/matrixA.m", "A");
        IMLOperatorToMatrix(projection_precond_operator, dense_precond_matrix);
        IMLVector iml_prec_b = projection_preconditioner.solve(iml_b);

        iml_x = dense_precond_matrix.solve_gmres(iml_prec_b);
//        std::cout << "b : \n";
//        std::cout << iml_b << std::endl;
//        std::cout << "precond b : \n";
//        std::cout << iml_prec_b << std::endl;
//        std::cout << "normal matrix " << dense_matrix;
//        std::cout << "precond matrix " << dense_precond_matrix;
//        std::cout << "x: \n";
//        std::cout << iml_x << std::endl;
//        std::cout << "check: \n";
//        std::cout << dense_precond_matrix * iml_x << std::endl;

        //iml_x = projection_preconditioner.solve(iml_b); int converged = 0;
        //std::cout << "error " << norm(projection_precond_operator * iml_b - iml_b) << std::endl;

//        int converged = GMRES(iml_operator, iml_x, iml_b, empty_preconditioner, iml_matrix, gmres_m, max_iter, tol);
//        int converged = GMRES(iml_operator, iml_x, iml_b, diag_preconditioner, iml_matrix, gmres_m, max_iter, tol);
//        int converged = GMRES(iml_operator, iml_x, iml_b, projection_preconditioner, iml_matrix, gmres_m, max_iter, tol);
//        int converged = GMRES(test_matrix, iml_x, projection_preconditioner.solve(iml_b), empty_preconditioner, iml_matrix, gmres_m, max_iter, tol);

//          int converged = GMRES(test_matrix, iml_x, diag_preconditioner.solve(iml_b), empty_preconditioner, iml_matrix, gmres_m, max_iter, tol);
//        int converged = GMRES(diag_precond_operator, iml_x, diag_preconditioner.solve(iml_b), empty_preconditioner, iml_matrix, gmres_m, max_iter, tol);
//        int converged = GMRES(projection_precond_operator, iml_x, projection_preconditioner.solve(iml_b), empty_preconditioner, iml_matrix, gmres_m, max_iter, tol);

//        int converged = CG(iml_operator, iml_x, iml_b, empty_preconditioner, max_iter, tol);
//        int converged = CG(iml_operator, iml_x, iml_b, diag_preconditioner, max_iter, tol);
//        int converged = CG(iml_operator, iml_x, iml_b, projection_preconditioner, max_iter, tol);

//        int converged = CG(projection_precond_operator, iml_x, projection_preconditioner.solve(iml_b), empty_preconditioner, max_iter, tol);
//        int converged = CG(diag_precond_operator, iml_x, diag_preconditioner.solve(iml_b), empty_preconditioner, max_iter, tol);
        iml_x.copy(&x);
        Hermes::Mixins::Loggable::Static::info("CG converged %d, tol %g, steps %d, ndofs %d", converged, tol, max_iter, ndof_ref);
        Solution<double>::vector_to_solution(x.v, ref_space, ref_sln);
    }
    cpu_time.tick();
    Hermes::Mixins::Loggable::Static::info("Solution: %g s", cpu_time.last());

    // Project the fine mesh solution onto the coarse mesh.
    Hermes::Mixins::Loggable::Static::info("Calculating error estimate and exact error.");
    OGProjection<double> ogProjection; ogProjection.project_global(space, ref_sln, sln);

    // Calculate element errors and total error estimate.
    adaptivity.set_space(space);
    errorCalculator.calculate_errors(sln, exact_sln, false);
    double err_coarse_exact_rel = errorCalculator.get_total_error_squared() * 100;

    errorCalculator.calculate_errors(ref_sln, exact_sln, false);
    double err_fine_exact_rel = errorCalculator.get_total_error_squared() * 100;

    errorCalculator.calculate_errors(sln, ref_sln, true);
    double err_est_rel = errorCalculator.get_total_error_squared() * 100;

    cpu_time.tick();
    Hermes::Mixins::Loggable::Static::info("Error calculation: %g s", cpu_time.last());
    
    // Report results.
    Hermes::Mixins::Loggable::Static::info("ndof_coarse: %d, ndof_fine: %d", space->get_num_dofs(), ref_space->get_num_dofs());
    Hermes::Mixins::Loggable::Static::info("err_est_rel: %g%%, err_coarse_exact_rel: %g%%, err_fine_exact_rel: %g%%", err_est_rel, err_coarse_exact_rel, err_fine_exact_rel);

    // Time measurement.
    cpu_time.tick();
    double accum_time = cpu_time.accumulated();
    
    // View the coarse mesh solution and polynomial orders.
//    oview.show(space);
    sview_coarse.show(sln);
    sview_fine.show(ref_sln);
    oview.show(ref_space);
    //getchar();

    // Add entry to DOF and CPU convergence graphs.
    graph_dof_est.add_values(space->get_num_dofs(), err_est_rel);
    graph_dof_est.save("conv_dof_est.dat");
    graph_cpu_est.add_values(accum_time, err_est_rel);
    graph_cpu_est.save("conv_cpu_est.dat");
    graph_dof_exact.add_values(space->get_num_dofs(), err_coarse_exact_rel);
    graph_dof_exact.save("conv_dof_exact.dat");
    graph_cpu_exact.add_values(accum_time, err_coarse_exact_rel);
    graph_cpu_exact.save("conv_cpu_exact.dat");
    
    cpu_time.tick(Hermes::Mixins::TimeMeasurable::HERMES_SKIP);

    // If err_est too large, adapt the mesh. The NDOF test must be here, so that the solution may be visualized
    // after ending due to this criterion.
    if (err_coarse_exact_rel < ERR_STOP)
      done = true;
    else
      done = adaptivity.adapt(&selector);
   
    cpu_time.tick();
    Hermes::Mixins::Loggable::Static::info("Adaptation: %g s", cpu_time.last());
    
    // Increase the counter of adaptivity steps.
    if (done == false)  
      as++;
  }
  while (done == false);

  Hermes::Mixins::Loggable::Static::info("Total running time: %g s", cpu_time.accumulated());

  // Wait for all views to be closed.
  Views::View::wait();
  return 0;
}
