#ifndef MYCG_H
#define MYCG_H

#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;

class CGSolver
{
public:
    int solve(SimpleVector<double>* x, SimpleVector<double>* b, double &tol, int &max_iter);
    virtual void apply_matrix(SimpleVector<double>* vector, SimpleVector<double>* result) = 0;
    virtual void solve_for_precond(SimpleVector<double>* vector, SimpleVector<double>* result) = 0;

    // the rest should not be here, just for visualisation
    Views::ScalarView* iterative_solution2;
    bool show_progress;

    SpaceSharedPtr<double> coarse_space;
    SpaceSharedPtr<double> fine_space;

};

#endif // MYCG_H
