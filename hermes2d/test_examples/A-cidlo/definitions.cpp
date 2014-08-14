#include "definitions.h"

/* Weak forms */

CustomWeakFormPoisson::CustomWeakFormPoisson(ProblemDefinition* definition, Perms perms, bool external_dirichlet_lift) : Hermes::Hermes2D::WeakForm<double>(1)
{
  // Matrix forms.
    add_matrix_form(new Hermes::Hermes2D::WeakFormsH1::DefaultMatrixFormDiffusion<double>(0, 0, definition->labels_air, new Hermes::Hermes1DFunction<double>(perms.EPS_AIR)));
    add_matrix_form(new Hermes::Hermes2D::WeakFormsH1::DefaultMatrixFormDiffusion<double>(0, 0, definition->labels_kartit, new Hermes::Hermes1DFunction<double>(perms.EPS_KARTIT)));
    add_matrix_form(new Hermes::Hermes2D::WeakFormsH1::DefaultMatrixFormDiffusion<double>(0, 0, definition->labels_full, new Hermes::Hermes1DFunction<double>(perms.EPS_FULL)));
    add_matrix_form(new Hermes::Hermes2D::WeakFormsH1::DefaultMatrixFormDiffusion<double>(0, 0, definition->labels_empty, new Hermes::Hermes1DFunction<double>(perms.EPS_EMPTY)));

    //Vector forms - for electrostatic should be removed(coeff = 0)
    if(definition->SOURCE_TERM != 0.0)
        add_vector_form(new Hermes::Hermes2D::WeakFormsH1::DefaultVectorFormVol<double>(0, Hermes::HERMES_ANY, new Hermes::Hermes2DFunction<double>(definition->SOURCE_TERM)));

    // Dirichlet lift is not done using standard Hermes tools.
    // It is calculated separately (using relative permitivity == 1 in whole domain) and than added to the final solution
    if(external_dirichlet_lift)
    {
        add_vector_form(new GradDirichletLiftTimesGradTest(0, definition->labels_air, -perms.EPS_AIR, 0));
        add_vector_form(new GradDirichletLiftTimesGradTest(0, definition->labels_kartit, -perms.EPS_KARTIT, 0));
        add_vector_form(new GradDirichletLiftTimesGradTest(0, definition->labels_full, -perms.EPS_FULL, 0));
        add_vector_form(new GradDirichletLiftTimesGradTest(0, definition->labels_empty, -perms.EPS_EMPTY, 0));
    }

}

CustomWeakFormPermitivity::CustomWeakFormPermitivity(ProblemDefinition* definition, Perms perms) : Hermes::Hermes2D::WeakForm<double>(1)
{
    // Matrix forms.
    add_matrix_form(new Hermes::Hermes2D::WeakFormsH1::DefaultMatrixFormVol<double>(0, 0));

    //Vector forms
    add_vector_form(new Hermes::Hermes2D::WeakFormsH1::DefaultVectorFormVol<double>(0, definition->labels_air, new Hermes::Hermes2DFunction<double>(perms.EPS_AIR/EPS0)));
    add_vector_form(new Hermes::Hermes2D::WeakFormsH1::DefaultVectorFormVol<double>(0, definition->labels_kartit, new Hermes::Hermes2DFunction<double>(perms.EPS_KARTIT/EPS0)));
    add_vector_form(new Hermes::Hermes2D::WeakFormsH1::DefaultVectorFormVol<double>(0, definition->labels_full, new Hermes::Hermes2DFunction<double>(perms.EPS_FULL/EPS0)));
    add_vector_form(new Hermes::Hermes2D::WeakFormsH1::DefaultVectorFormVol<double>(0, definition->labels_empty, new Hermes::Hermes2DFunction<double>(perms.EPS_EMPTY/EPS0)));

}


WeakFormChangingPermInFull::WeakFormChangingPermInFull(const PGDSolutions* pgd_sols) : Hermes::Hermes2D::WeakForm<double>(1)
{
    assert(pgd_sols->parameters[0].size() == pgd_sols->solutions.size());

    double w1_air = pgd_sols->perms.EPS_AIR * pgd_sols->actual_parameter[0].int_F_F();
    double w1_kartit = pgd_sols->perms.EPS_KARTIT * pgd_sols->actual_parameter[0].int_F_F();
    double w1_empty = pgd_sols->perms.EPS_EMPTY * pgd_sols->actual_parameter[0].int_F_F();

    // integral with x -- full
    double w1_full = pgd_sols->actual_parameter[0].int_x_F_F();

    add_matrix_form(new WeakFormsH1::DefaultMatrixFormDiffusion<double>(0, 0, pgd_sols->definition->labels_air, new Hermes1DFunction<double>(w1_air)));
    add_matrix_form(new WeakFormsH1::DefaultMatrixFormDiffusion<double>(0, 0, pgd_sols->definition->labels_kartit, new Hermes1DFunction<double>(w1_kartit)));
    add_matrix_form(new WeakFormsH1::DefaultMatrixFormDiffusion<double>(0, 0, pgd_sols->definition->labels_full, new Hermes1DFunction<double>(w1_full)));
    add_matrix_form(new WeakFormsH1::DefaultMatrixFormDiffusion<double>(0, 0, pgd_sols->definition->labels_empty, new Hermes1DFunction<double>(w1_empty)));

    // Residual forms.
    std::vector<double> w2_air, w2_kartit, w2_empty, w2_full;

    std::cout << "pushing " << pgd_sols->solutions.size() << " coefficients" << std::endl;
    for(int i = 0; i < pgd_sols->solutions.size(); i++)
    {
        w2_air.push_back(-pgd_sols->perms.EPS_AIR * pgd_sols->actual_parameter[0].int_F_ExtF(pgd_sols->parameters[0][i]));
        w2_kartit.push_back(-pgd_sols->perms.EPS_KARTIT * pgd_sols->actual_parameter[0].int_F_ExtF(pgd_sols->parameters[0][i]));
        w2_empty.push_back(-pgd_sols->perms.EPS_EMPTY * pgd_sols->actual_parameter[0].int_F_ExtF(pgd_sols->parameters[0][i]));

        // integral with x -- full
        w2_full.push_back(-pgd_sols->actual_parameter[0].int_x_F_ExtF(pgd_sols->parameters[0][i]));
    }

    add_vector_form(new GradPreviousSolsTimesGradTest(0, pgd_sols->definition->labels_air, w2_air));
    add_vector_form(new GradPreviousSolsTimesGradTest(0, pgd_sols->definition->labels_kartit, w2_kartit));
    add_vector_form(new GradPreviousSolsTimesGradTest(0, pgd_sols->definition->labels_full, w2_full));
    add_vector_form(new GradPreviousSolsTimesGradTest(0, pgd_sols->definition->labels_empty, w2_empty));

    // RHS residual forms
    if(pgd_sols->definition->SOURCE_TERM != 0.0)
    {
        double w3 = pgd_sols->actual_parameter[0].int_F();
        add_vector_form(new WeakFormsH1::DefaultVectorFormVol<double>(0, Hermes::HERMES_ANY, new Hermes::Hermes2DFunction<double>(w3 * pgd_sols->definition->SOURCE_TERM)));
    }

    // Dirichlet lift
    if(pgd_sols->definition->use_dirichlet_lift())
    {
        double w3_air = -pgd_sols->perms.EPS_AIR * pgd_sols->actual_parameter[0].int_F();
        double w3_empty = -pgd_sols->perms.EPS_EMPTY * pgd_sols->actual_parameter[0].int_F();
        double w3_kartit = -pgd_sols->perms.EPS_KARTIT * pgd_sols->actual_parameter[0].int_F();

        double w3_full = -pgd_sols->actual_parameter[0].int_x_F();

        int dirichlet_lift_idx = pgd_sols->solutions.size();
        add_vector_form(new GradDirichletLiftTimesGradTest(0, pgd_sols->definition->labels_air, w3_air, dirichlet_lift_idx));
        add_vector_form(new GradDirichletLiftTimesGradTest(0, pgd_sols->definition->labels_empty, w3_empty, dirichlet_lift_idx));
        add_vector_form(new GradDirichletLiftTimesGradTest(0, pgd_sols->definition->labels_kartit, w3_kartit, dirichlet_lift_idx));
        add_vector_form(new GradDirichletLiftTimesGradTest(0, pgd_sols->definition->labels_full, w3_full, dirichlet_lift_idx));
    }

    Hermes::vector<MeshFunctionSharedPtr<double> > previous_sols;
    for(int i = 0; i < pgd_sols->solutions.size(); i++)
        previous_sols.push_back(pgd_sols->solutions.at(i));
    if(pgd_sols->definition->use_dirichlet_lift())
    {
        previous_sols.push_back(pgd_sols->dirichlet_lift);
        assert(previous_sols.size() == pgd_sols->solutions.size() + 1);
    }
    set_ext(previous_sols);
}

WeakFormMultipleColumns::WeakFormMultipleColumns(const PGDSolutions *pgd_sols)
{
    ProblemDefinition* problemDefinition = pgd_sols->definition;

    // zatim jen pro jeden celkovy sloupec
    assert(N_WIDTH_COARSE == 1);
    for(int column_idx = 0; column_idx < N_WIDTH_COARSE; column_idx++)
    {
        for(int row_idx = 0; row_idx < N_HEIGHT_COARSE; row_idx++)
        {
            Hermes::vector<std::string> element_markers = problemDefinition->labels_square[column_idx][row_idx];
            double w1 = pgd_sols->actual_parameter[0].int_epsx_F_F(row_idx, pgd_sols->perms);
            add_matrix_form(new WeakFormsH1::DefaultMatrixFormDiffusion<double>(0, 0, element_markers, new Hermes1DFunction<double>(w1)));

        }
    }

}

GradPreviousSolsTimesGradTest::GradPreviousSolsTimesGradTest(int i, Hermes::vector<std::string> areas, std::vector<double> coeffs)
    : VectorFormVol<double>(i), coeffs(coeffs)
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
        for(int prev_idx = 0; prev_idx < coeffs.size(); prev_idx++)
        {
            result += wt[i] * coeffs[prev_idx] * (ext[prev_idx]->dx[i] * v->dx[i] + ext[prev_idx]->dy[i] * v->dy[i]);
        }
    }
    return result;
}

Ord GradPreviousSolsTimesGradTest::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                                      Geom<Ord> *e, Func<Ord> **ext) const
{
    Ord result = Ord(0);
    for (int i = 0; i < n; i++) {
        for(int prev_idx = 0; prev_idx < coeffs.size(); prev_idx++)
        {
            result += wt[i] * coeffs[prev_idx] * (ext[prev_idx]->dx[i] * v->dx[i] + ext[prev_idx]->dy[i] * v->dy[i]);
        }
    }
    return result;
}

VectorFormVol<double>* GradPreviousSolsTimesGradTest::clone() const
{
    return new GradPreviousSolsTimesGradTest(this->i, this->areas, this->coeffs);
}

GradDirichletLiftTimesGradTest::GradDirichletLiftTimesGradTest(int i, Hermes::vector<std::string> areas, double coeff, int extIdx)
    : VectorFormVol<double>(i), ext_idx(extIdx), coeff(coeff)
{
    this->set_areas(areas);
}

GradDirichletLiftTimesGradTest::~GradDirichletLiftTimesGradTest()
{
}

double GradDirichletLiftTimesGradTest::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v,
                                           Geom<double> *e, Func<double> **ext) const
{
    double result = 0;
    for (int i = 0; i < n; i++) {
        result += wt[i] * (ext[ext_idx]->dx[i] * v->dx[i] + ext[ext_idx]->dy[i] * v->dy[i]);
    }
    return coeff * result;
}

Ord GradDirichletLiftTimesGradTest::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                                      Geom<Ord> *e, Func<Ord> **ext) const
{
    Ord result = Ord(0);
    for (int i = 0; i < n; i++) {
        result += wt[i] * (ext[ext_idx]->dx[i] * v->dx[i] + ext[ext_idx]->dy[i] * v->dy[i]);
    }
    return coeff * result;
}

VectorFormVol<double>* GradDirichletLiftTimesGradTest::clone() const
{
    return new GradDirichletLiftTimesGradTest(this->i, this->areas, this->coeff, this->ext_idx);
}

