#include "hermes2d.h"
#include "problem.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;

class CustomWeakFormPoisson : public Hermes::Hermes2D::WeakForm<double>
{
public:
    CustomWeakFormPoisson(ProblemDefinition definition, Perms perms);
};

class CustomWeakFormPermitivity : public Hermes::Hermes2D::WeakForm<double>
{
public:
    CustomWeakFormPermitivity(ProblemDefinition definition, Perms perms);
};

struct PGDSolutions
{
    std::vector<Function1D> parameters;
    std::vector<MeshFunctionSharedPtr<double> > solutions;

    double get_pt_value(double x, double y, double eps)
    {
        assert(solutions.size() == parameters.size());
        double result = 0;
        for(int i = 0; i < solutions.size(); i++)
        {
            result += parameters.at(i).value(eps) * solutions.at(i)->get_pt_value(x, y)->val[0];
        }

        return result;
    }
};

class WeakFormChangingPermInFull : public Hermes::Hermes2D::WeakForm<double>
{
public:
    WeakFormChangingPermInFull(ProblemDefinition definition, Perms perms, PGDSolutions pgd_solutions);
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
  int idx_i;
  std::vector<double> coeffs;
};

class EnergyIntegralCalculator : public PostProcessing::VolumetricIntegralCalculator<double>
{
public:
    EnergyIntegralCalculator(MeshFunctionSharedPtr<double> source_function, ProblemDefinition definition, Perms perms) :
        PostProcessing::VolumetricIntegralCalculator<double>(source_function, 1), definition(definition)
    {
        //memset(label_to_eps, 0, MAX_LABELS * sizeof(double));
        for(int i = 0; i < MAX_LABELS; i++)
            label_to_eps[i] = 0.;

        add_labels(definition.labels_air, perms.EPS_AIR);
        add_labels(definition.labels_kartit, perms.EPS_KARTIT);
        add_labels(definition.labels_full, perms.EPS_FULL);
        add_labels(definition.labels_empty, perms.EPS_EMPTY);
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

    ProblemDefinition definition;
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
