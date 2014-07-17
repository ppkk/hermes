#include "definitions.h"

/* Weak forms */

CustomWeakFormPoisson::CustomWeakFormPoisson(ProblemDefinition definition, Perms perms) : Hermes::Hermes2D::WeakForm<double>(1)
{
  // Jacobian forms.
    add_matrix_form(new Hermes::Hermes2D::WeakFormsH1::DefaultMatrixFormDiffusion<double>(0, 0, definition.labels_air, new Hermes::Hermes1DFunction<double>(perms.EPS_AIR)));
    add_matrix_form(new Hermes::Hermes2D::WeakFormsH1::DefaultMatrixFormDiffusion<double>(0, 0, definition.labels_kartit, new Hermes::Hermes1DFunction<double>(perms.EPS_KARTIT)));
    add_matrix_form(new Hermes::Hermes2D::WeakFormsH1::DefaultMatrixFormDiffusion<double>(0, 0, definition.labels_full, new Hermes::Hermes1DFunction<double>(perms.EPS_FULL)));
    add_matrix_form(new Hermes::Hermes2D::WeakFormsH1::DefaultMatrixFormDiffusion<double>(0, 0, definition.labels_empty, new Hermes::Hermes1DFunction<double>(perms.EPS_EMPTY)));
}

CustomWeakFormPermitivity::CustomWeakFormPermitivity(ProblemDefinition definition, Perms perms) : Hermes::Hermes2D::WeakForm<double>(1)
{
  // Jacobian forms.
    add_matrix_form(new Hermes::Hermes2D::WeakFormsH1::DefaultMatrixFormVol<double>(0, 0));
    add_vector_form(new Hermes::Hermes2D::WeakFormsH1::DefaultVectorFormVol<double>(0, definition.labels_air, new Hermes::Hermes2DFunction<double>(perms.EPS_AIR/EPS0)));
    add_vector_form(new Hermes::Hermes2D::WeakFormsH1::DefaultVectorFormVol<double>(0, definition.labels_kartit, new Hermes::Hermes2DFunction<double>(perms.EPS_KARTIT/EPS0)));
    add_vector_form(new Hermes::Hermes2D::WeakFormsH1::DefaultVectorFormVol<double>(0, definition.labels_full, new Hermes::Hermes2DFunction<double>(perms.EPS_FULL/EPS0)));
    add_vector_form(new Hermes::Hermes2D::WeakFormsH1::DefaultVectorFormVol<double>(0, definition.labels_empty, new Hermes::Hermes2DFunction<double>(perms.EPS_EMPTY/EPS0)));

}


WeakFormChangingPermInFull::WeakFormChangingPermInFull(ProblemDefinition definition, Perms perms, PGDSolutions pgd_solutions) : Hermes::Hermes2D::WeakForm<double>(1)
{
    assert(pgd_solutions.parameters.size() == pgd_solutions.solutions.size() + 1);

    // Jacobian forms.
    double w1_air = perms.EPS_AIR * pgd_solutions.parameters.back().int_F_F();
    double w1_kartit = perms.EPS_KARTIT * pgd_solutions.parameters.back().int_F_F();
    double w1_empty = perms.EPS_EMPTY * pgd_solutions.parameters.back().int_F_F();
    double w1_full = pgd_solutions.parameters.back().int_x_F_F();
    add_matrix_form(new WeakFormsH1::DefaultMatrixFormDiffusion<double>(0, 0, definition.labels_air, new Hermes1DFunction<double>(w1_air)));
    add_matrix_form(new WeakFormsH1::DefaultMatrixFormDiffusion<double>(0, 0, definition.labels_kartit, new Hermes1DFunction<double>(w1_kartit)));
    add_matrix_form(new WeakFormsH1::DefaultMatrixFormDiffusion<double>(0, 0, definition.labels_full, new Hermes1DFunction<double>(w1_full)));
    add_matrix_form(new WeakFormsH1::DefaultMatrixFormDiffusion<double>(0, 0, definition.labels_empty, new Hermes1DFunction<double>(w1_empty)));

    // Residual forms.
    std::vector<double> w2_air;
    std::vector<double> w2_kartit;
    std::vector<double> w2_empty;
    std::vector<double> w2_full;

    for(int i = 0; i < pgd_solutions.solutions.size(); i++)
    {
        w2_air.push_back(perms.EPS_AIR * pgd_solutions.parameters.back().int_F_ExtF(pgd_solutions.parameters[i]));
        w2_kartit.push_back(perms.EPS_KARTIT * pgd_solutions.parameters.back().int_F_ExtF(pgd_solutions.parameters[i]));
        w2_empty.push_back(perms.EPS_EMPTY * pgd_solutions.parameters.back().int_F_ExtF(pgd_solutions.parameters[i]));
        w2_full.push_back(pgd_solutions.parameters.back().int_x_F_ExtF(pgd_solutions.parameters[i]));
    }

    add_vector_form(new GradPreviousSolsTimesGradTest(0, definition.labels_air, w2_air));
    add_vector_form(new GradPreviousSolsTimesGradTest(0, definition.labels_kartit, w2_kartit));
    add_vector_form(new GradPreviousSolsTimesGradTest(0, definition.labels_full, w2_full));
    add_vector_form(new GradPreviousSolsTimesGradTest(0, definition.labels_empty, w2_empty));
}

GradPreviousSolsTimesGradTest::GradPreviousSolsTimesGradTest(int i, Hermes::vector<std::string> areas, std::vector<double> coeffs)
    : VectorFormVol<double>(i), idx_i(i), coeffs(coeffs)
{
    this->set_areas(areas);
}

GradPreviousSolsTimesGradTest::~GradPreviousSolsTimesGradTest()
{
}

double GradPreviousSolsTimesGradTest::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v,
                                           Geom<double> *e, Func<double> **ext) const
{
    double result = 0;
    for (int i = 0; i < n; i++) {
        double value = 0;
        for(int prev_idx = 0; prev_idx < coeffs.size(); prev_idx++)
        {
            value += coeffs[prev_idx] * (ext[prev_idx]->dx[i] * v->dx[i] + ext[prev_idx]->dy[i] * v->dy[i]);
        }
        result += wt[i] * value;
    }
    //   Here is the - sign
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!
    return - result;
}

Ord GradPreviousSolsTimesGradTest::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                                      Geom<Ord> *e, Func<Ord> **ext) const
{
    Ord result = Ord(0);
    for (int i = 0; i < n; i++) {
        Ord value = Ord(0);
        for(int prev_idx = 0; prev_idx < coeffs.size(); prev_idx++)
        {
            value += coeffs[prev_idx] * (ext[prev_idx]->dx[i] * v->dx[i] + ext[prev_idx]->dy[i] * v->dy[i]);
        }
        result += wt[i] * value;
    }
    //   Here is the - sign
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!
    return - result;
}

VectorFormVol<double>* GradPreviousSolsTimesGradTest::clone() const
{
    return new GradPreviousSolsTimesGradTest(this->i, this->areas, this->coeffs);
}

