// This file is part of HermesCommon
//
// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Email: hpfem-group@unr.edu, home page: http://hpfem.org/.
//
// Hermes2D is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation; either version 2 of the License,
// or (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
/*! \file nonlinear_solver.h
\brief General nonlinear solver functionality.
*/
#ifndef __HERMES_COMMON_LINEAR_SOLVER_H_
#define __HERMES_COMMON_LINEAR_SOLVER_H_

#include "discrete_problem_linear.h"

namespace Hermes
{
  namespace Hermes2D
  {
    /// \brief Base class for defining interface for nonlinear solvers.
    ///
    template <typename Scalar>
    class LinearSolver : public Hermes::Mixins::Loggable
    {
    public:
      LinearSolver(DiscreteProblemLinear<Scalar>* dp);

      ~LinearSolver();

      /// Basic solve method.
      virtual void solve();

      Scalar *get_sln_vector();
    protected:
      DiscreteProblemLinear<Scalar>* dp; ///< FE problem being solved.

      /// The solution vector.
      Scalar* sln_vector;

      /// Jacobian.
      SparseMatrix<Scalar>* jacobian;

      /// Residual.
      Vector<Scalar>* residual;

      /// Linear solver.
      LinearMatrixSolver<Scalar>* matrix_solver;
    };
  }
}
#endif