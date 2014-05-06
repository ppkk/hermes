#include "mycg.h"

int CGSolver::solve(SimpleVector<double>* x, SimpleVector<double>* b, double &tol, int &max_iter)
{
    assert(b->get_size() == x->get_size());
    double beta = 1, alpha = 1, rho = 1, rho_1 = 1;
    SimpleVector<double> r, Ax, z, p, q;
    r.alloc(b->get_size());
    z.alloc(b->get_size());
    p.alloc(b->get_size());
    q.alloc(b->get_size());
    Ax.alloc(b->get_size());
    r.set_vector(b);
    apply_matrix(x, &Ax);
    r.add_vector_multiple(&Ax, -1);

    double norm_b = get_l2_norm(b);
    if(norm_b == 0.0)
        norm_b = 1.0;

    double resid = get_l2_norm(&r) / norm_b;

    if (resid <= tol) {
       tol = resid;
       max_iter = 0;
       return 0;
     }

    for (int i = 1; i <= max_iter; i++) {

        if(show_progress)
        {
            MeshFunctionSharedPtr<double> tmp_solution(new Solution<double>());
            Solution<double>::vector_to_solution(x->v, fine_space, tmp_solution);


            double minr = -2, maxr = 2;
            Views::LinearizerCriterionFixed lincrit(4);iterative_solution2 = new Views::ScalarView("Iterative solution", new Views::WinGeom(0, 0, 440, 350));
            iterative_solution2->show_mesh(true);
            iterative_solution2->fix_scale_width(50);
            iterative_solution2->set_min_max_range(minr,maxr);
            iterative_solution2->get_linearizer()->set_criterion(lincrit);

            iterative_solution2->show(tmp_solution);
        }


        solve_for_precond(&r, &z);
        rho = get_dot_product(&r, &z);

        if (i == 1)
          p.set_vector(&z);
        else {
          beta = rho / rho_1;
          p.set_vector(&z);
          p.add_vector_multiple(&p, beta);
        }

        apply_matrix(&p, &q);
        alpha = rho / get_dot_product(&p, &q);

        x->add_vector_multiple(&p, alpha);
        r.add_vector_multiple(&q, -alpha);

        resid = get_l2_norm(&r) / norm_b;
        std::cout << "residual " << resid << std::endl;
        if (resid <= tol) {
          tol = resid;
          max_iter = i;
          return 0;
        }

        rho_1 = rho;
      }

      tol = resid;
      return 1;
}
