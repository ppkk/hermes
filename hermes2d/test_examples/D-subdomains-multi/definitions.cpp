#include "definitions.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::WeakFormsH1;

CustomWeakFormCurrentHeat::CustomWeakFormCurrentHeat() : Hermes::Hermes2D::WeakForm<double>(2)
{
  this->add_matrix_form(new DefaultJacobianDiffusion<double>(0, 0));
  this->add_matrix_form(new DefaultJacobianDiffusion<double>(1, 1));

  // Residual.
  this->add_vector_form(new DefaultResidualDiffusion<double>(0));
  this->add_vector_form(new DefaultVectorFormVol<double>(0));
  this->add_vector_form(new DefaultResidualDiffusion<double>(1));
  this->add_vector_form(new DefaultVectorFormVol<double>(1));
//  this->add_vector_form(new CustomVectorFormVol(2));
//  this->add_vector_form(new DefaultResidualDiffusion<double>(2));
//  this->add_vector_form(new DefaultVectorFormVol<double>(2));
}

CustomWeakFormCurrent::CustomWeakFormCurrent() : Hermes::Hermes2D::WeakForm<double>(1)
{
  this->add_matrix_form(new DefaultJacobianDiffusion<double>(0, 0));

  // Residual.
  this->add_vector_form(new DefaultResidualDiffusion<double>(0));
  this->add_vector_form(new DefaultVectorFormVol<double>(0));
}

CustomWeakFormHeat::CustomWeakFormHeat(Hermes::Hermes2D::Solution<double>* solution_current) : Hermes::Hermes2D::WeakForm<double>(1)
{
  this->add_matrix_form(new DefaultJacobianDiffusion<double>(0, 0));

  // Residual.
  this->add_vector_form(new DefaultResidualDiffusion<double>(0));
  this->add_vector_form(new DefaultVectorFormVol<double>(0));
  CustomVectorFormVol *custom = new CustomVectorFormVol(0, solution_current);
  this->add_vector_form(custom);
}

CustomVectorFormVol::CustomVectorFormVol(int i, Hermes::Hermes2D::Solution<double>* solution_current)
  : VectorFormVol<double>(i)
{
    m_solution_current = solution_current;
    ext.push_back(solution_current);
}

double CustomVectorFormVol::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v,
  Geom<double> *e, ExtData<double> *ext) const
{
  double result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * ext->fn[0]->val[i] * v->val[i];
  return result;
}

Ord CustomVectorFormVol::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
  Geom<Ord> *e, ExtData<Ord> *ext) const
{
  Ord result = Ord(0);
  for (int i = 0; i < n; i++)
    result += wt[i] * ext->fn[0]->val[i] * v->val[i];
  return result;
}

VectorFormVol<double>* CustomVectorFormVol::clone()
{
    return new CustomVectorFormVol(1, m_solution_current);
}
