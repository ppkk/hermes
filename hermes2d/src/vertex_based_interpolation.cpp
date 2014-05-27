// This file is part of Hermes2D.
//
// Hermes2D is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D.  If not, see <http://www.gnu.org/licenses/>.

#include "vertex_based_interpolation.h"
#include "forms.h"

namespace Hermes
{
  namespace Hermes2D
  {
    template<typename Scalar>
    void VertexBasedInterpolation<Scalar>::interpolate(SpaceSharedPtr<Scalar> src_space, Scalar* src_sln_vector, SpaceSharedPtr<Scalar> out_space, Scalar*& out_sln_vector, bool include_dirichlet)
    {
      out_sln_vector = calloc_with_check<Scalar>(out_space->get_num_dofs());
      MeshFunctionSharedPtr<Scalar> src_fn(new Solution<Scalar>(src_space->get_mesh()));
      Solution<Scalar>::vector_to_solution(src_sln_vector, src_space, src_fn, include_dirichlet);

      // Initialize states && previous iterations.
      int num_states;
      Traverse::State** states;
      Hermes::vector<MeshSharedPtr> meshes(src_space->get_mesh(), out_space->get_mesh());
      Traverse trav(2);
      states = trav.get_states(meshes, num_states);

      for (int state_i = 0; state_i < num_states; state_i++)
      {
        Traverse::State* current_state = states[state_i];

        for (unsigned char edge = 0; edge < current_state->rep->nvert; edge++)
        {
          typename Space<Scalar>::NodeData* out_nd = out_space->ndata + current_state->e[1]->vn[edge]->id;
          if (out_nd->vertex_bc_coef != nullptr)
            continue;
          if (current_state->e[0]->vn[edge]->id == current_state->e[1]->vn[edge]->id && out_space->ndata[current_state->e[1]->vn[edge]->id].n == 1 && src_space->ndata[current_state->e[0]->vn[edge]->id].n == 1)
          {
            typename Space<Scalar>::NodeData* src_nd = src_space->ndata + current_state->e[0]->vn[edge]->id;
            if (src_space->ndata[current_state->e[0]->vn[edge]->id].n == 1 && src_nd->dof >= 0)
              out_sln_vector[out_nd->dof] = src_sln_vector[src_nd->dof];
          }
          else
          {
            Func<Scalar>* func = (dynamic_cast<Solution<Scalar>*>(src_fn.get()))->get_pt_value(current_state->e[1]->vn[edge]->x, current_state->e[1]->vn[edge]->y);
            out_sln_vector[out_nd->dof] = func->val[0];
          }
        }
      }
    }

    template class HERMES_API VertexBasedInterpolation<double>;
    template class HERMES_API VertexBasedInterpolation<std::complex<double> >;
  }
}
