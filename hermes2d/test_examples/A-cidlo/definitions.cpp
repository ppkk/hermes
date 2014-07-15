#include "definitions.h"

/* Weak forms */

CustomWeakFormPoisson::CustomWeakFormPoisson(ProblemDefinition definition) : Hermes::Hermes2D::WeakForm<double>(1)
{
  // Jacobian forms.
    add_matrix_form(new Hermes::Hermes2D::WeakFormsH1::DefaultMatrixFormDiffusion<double>(0, 0, definition.labels_air, new Hermes::Hermes1DFunction<double>(definition.EPS_AIR)));
    add_matrix_form(new Hermes::Hermes2D::WeakFormsH1::DefaultMatrixFormDiffusion<double>(0, 0, definition.labels_kartit, new Hermes::Hermes1DFunction<double>(definition.EPS_KARTIT)));
    add_matrix_form(new Hermes::Hermes2D::WeakFormsH1::DefaultMatrixFormDiffusion<double>(0, 0, definition.labels_full, new Hermes::Hermes1DFunction<double>(definition.EPS_FULL)));
    add_matrix_form(new Hermes::Hermes2D::WeakFormsH1::DefaultMatrixFormDiffusion<double>(0, 0, definition.labels_empty, new Hermes::Hermes1DFunction<double>(definition.EPS_EMPTY)));
};

CustomWeakFormPermitivity::CustomWeakFormPermitivity(ProblemDefinition definition) : Hermes::Hermes2D::WeakForm<double>(1)
{
  // Jacobian forms.
    add_matrix_form(new Hermes::Hermes2D::WeakFormsH1::DefaultMatrixFormVol<double>(0, 0));
    add_vector_form(new Hermes::Hermes2D::WeakFormsH1::DefaultVectorFormVol<double>(0, definition.labels_air, new Hermes::Hermes2DFunction<double>(definition.EPS_AIR/EPS0)));
    add_vector_form(new Hermes::Hermes2D::WeakFormsH1::DefaultVectorFormVol<double>(0, definition.labels_kartit, new Hermes::Hermes2DFunction<double>(definition.EPS_KARTIT/EPS0)));
    add_vector_form(new Hermes::Hermes2D::WeakFormsH1::DefaultVectorFormVol<double>(0, definition.labels_full, new Hermes::Hermes2DFunction<double>(definition.EPS_FULL/EPS0)));
    add_vector_form(new Hermes::Hermes2D::WeakFormsH1::DefaultVectorFormVol<double>(0, definition.labels_empty, new Hermes::Hermes2DFunction<double>(definition.EPS_EMPTY/EPS0)));

};
