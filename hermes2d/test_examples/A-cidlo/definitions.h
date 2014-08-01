#include "hermes2d.h"
#include "discretization.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;

class CustomWeakFormPoisson : public Hermes::Hermes2D::WeakForm<double>
{
public:
    CustomWeakFormPoisson(ProblemDefinition* definition, Perms perms, bool external_dirichlet_lift);
};

class CustomWeakFormPermitivity : public Hermes::Hermes2D::WeakForm<double>
{
public:
    CustomWeakFormPermitivity(ProblemDefinition* definition, Perms perms);
};

class WeakFormChangingPermInFull : public Hermes::Hermes2D::WeakForm<double>
{
public:
    WeakFormChangingPermInFull(const PGDSolutions* pgd_sols);
};

class GradPreviousSolsTimesGradTest : public VectorFormVol<double>
{
public:
  GradPreviousSolsTimesGradTest(int i, Hermes::vector<std::string> areas, std::vector<double> coeffs);

  ~GradPreviousSolsTimesGradTest();

  virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v,
    Geom<double> *e, Func<double> **ext) const;

  virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *v,
    Geom<Hermes::Ord> *e, Func<Ord> **ext) const;

  virtual VectorFormVol<double>* clone() const;

private:
  std::vector<double> coeffs;
};

class GradDirichletLiftTimesGradTest : public VectorFormVol<double>
{
public:
  GradDirichletLiftTimesGradTest(int i, Hermes::vector<std::string> areas, double coeff, int ext_idx);

  ~GradDirichletLiftTimesGradTest();

  virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v,
    Geom<double> *e, Func<double> **ext) const;

  virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *v,
    Geom<Hermes::Ord> *e, Func<Ord> **ext) const;

  virtual VectorFormVol<double>* clone() const;

private:
  int ext_idx;
  double coeff;
};

class EnergyIntegralCalculator : public PostProcessing::VolumetricIntegralCalculator<double>
{
public:
    EnergyIntegralCalculator(MeshFunctionSharedPtr<double> source_function, ProblemDefinition* definition, Perms perms) :
        PostProcessing::VolumetricIntegralCalculator<double>(source_function, 1), definition(definition)
    {
        //memset(label_to_eps, 0, MAX_LABELS * sizeof(double));
        for(int i = 0; i < MAX_LABELS; i++)
            label_to_eps[i] = 0.;

        add_labels(definition->labels_air, perms.EPS_AIR);
        add_labels(definition->labels_kartit, perms.EPS_KARTIT);
        add_labels(definition->labels_full, perms.EPS_FULL);
        add_labels(definition->labels_empty, perms.EPS_EMPTY);
    }

    void add_labels(Hermes::vector<std::string> labels, double eps)
    {
        for(int i = 0; i < labels.size(); i++)
        {

            Hermes::Hermes2D::Mesh::MarkersConversion::IntValid iv = source_functions[0]->get_mesh()->get_element_markers_conversion().get_internal_marker(labels.at(i));
            assert(iv.valid);
            assert(iv.marker < MAX_LABELS);
            //std::cout << "adding " << iv.marker << std::endl;
            label_to_eps[iv.marker] = eps;
        }

    }

    virtual void integral(int n, double* wt, Func<double> **fns, Geom<double> *e, double* result)
    {
        assert(label_to_eps[e->elem_marker] != 0.0);
        //if(label_to_eps[e->elem_marker] == 0.0)
        //    std::cout << "integrating " << e->elem_marker << std::endl;
        for (int i = 0; i < n; i++){
            result[0] += wt[i] * 0.5* (fns[0]->dx[i]* fns[0]->dx[i] + fns[0]->dy[i]* fns[0]->dy[i]) * label_to_eps[e->elem_marker];
        }
    };

    virtual void order(Func<Hermes::Ord> **fns, Hermes::Ord* result) {
        result[0] = Hermes::Ord(21);
    }

    ProblemDefinition* definition;
    static const int MAX_LABELS = 5000;
    double label_to_eps[MAX_LABELS];
};

class GradUGradVIngegralCalculator : public PostProcessing::VolumetricIntegralCalculator<double>
{
public:
    GradUGradVIngegralCalculator(Hermes::vector<MeshFunctionSharedPtr<double> > source_functions) :
        PostProcessing::VolumetricIntegralCalculator<double>(source_functions, 1)
    {
    }

    virtual void integral(int n, double* wt, Func<double> **fns, Geom<double> *e, double* result)
    {
        for (int i = 0; i < n; i++){
            result[0] += wt[i] * (fns[0]->dx[i]* fns[1]->dx[i] + fns[0]->dy[i]* fns[1]->dy[i]);
        }
    };

    virtual void order(Func<Hermes::Ord> **fns, Hermes::Ord* result) {
        result[0] = Hermes::Ord(21);
    }

};

class U_times_f_IngegralCalculator : public PostProcessing::VolumetricIntegralCalculator<double>
{
public:
    U_times_f_IngegralCalculator(MeshFunctionSharedPtr<double> source_function, double coeff) :
        PostProcessing::VolumetricIntegralCalculator<double>(source_function, 1), coeff(coeff)
    {
    }

    virtual void integral(int n, double* wt, Func<double> **fns, Geom<double> *e, double* result)
    {
        for (int i = 0; i < n; i++){
            result[0] += wt[i] * (fns[0]->val[i]);
        }
        result[0] *= coeff;
    }

    virtual void order(Func<Hermes::Ord> **fns, Hermes::Ord* result) {
        result[0] = Hermes::Ord(21);
    }

    double coeff;

};

class CombinationFilter : public Hermes::Hermes2D::DXDYFilter<double>
{
public:
    CombinationFilter(Hermes::vector<MeshFunctionSharedPtr<double> > slns, std::vector<std::vector<Function1D> > parameters, double parameter_value) :
        DXDYFilter(slns), parameters(parameters), parameter_value(parameter_value) {}
    virtual void filter_fn(int n, double* x, double* y, Hermes::vector<const double *> values, Hermes::vector<const double *> dx, Hermes::vector<const double *> dy, double* rslt, double* rslt_dx, double* rslt_dy)
    {
        // number of different parameters
        int num_parameters = parameters.size();

        // number of modes - steps in iterative improvement of the solution
        int num_modes = parameters[0].size();

        for(int i = 0; i < num_parameters; i++)
            assert(parameters[i].size() == num_modes);

        double parameter_values_precalc[num_parameters][num_modes];
        for(int parameter = 0; parameter < num_parameters; parameter++)
            for(int mode = 0; mode < num_modes; mode++)
                parameter_values_precalc[parameter][mode] = parameters[parameter][mode].value(parameter_value);

        for (int i = 0; i < n; i++)
        {
            rslt[i] = 0;
            rslt_dx[i] = 0;
            rslt_dy[i] = 0;


            // the last value is the Dirichlet lift
            assert(values.size() == num_modes + 1);

            for (int mode = 0; mode < num_modes; mode++)
            {
                double val = values.at(mode)[i];
                double val_dx = dx.at(mode)[i];
                double val_dy = dy.at(mode)[i];

                for(int parameter = 0; parameter < num_parameters; parameter++)
                {
                    val *= parameter_values_precalc[parameter][mode];
                    val_dx *= parameter_values_precalc[parameter][mode];
                    val_dy *= parameter_values_precalc[parameter][mode];
                }

                rslt[i] += val;
                rslt_dx[i] += val_dx;
                rslt_dy[i] += val_dy;
            }

            rslt[i] += values.at(num_modes)[i];
            rslt_dx[i] +=  dx.at(num_modes)[i];
            rslt_dy[i] +=  dy.at(num_modes)[i];
        }
    }
    
    virtual MeshFunction<double>* clone() const
    {
        Hermes::vector<MeshFunctionSharedPtr<double> > slns;

        assert(this->num == parameters[0].size() + 1);
        for (int i = 0; i < this->num; i++)
        {
          slns.push_back(this->sln[i]->clone());
        }

        CombinationFilter* filter = new CombinationFilter(slns, parameters, parameter_value);
        return filter;
    }

private:
    std::vector<std::vector<Function1D> > parameters;
    double parameter_value;
};
