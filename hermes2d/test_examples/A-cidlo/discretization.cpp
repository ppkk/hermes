#include "discretization.h"
#include "definitions.h"

MeshFunctionSharedPtr<double> PGDSolutions::get_filter(double parameter_value, int num_modes)
{
    if(num_modes == -1)
        num_modes = solutions.size();
    assert((num_modes >= 0) && (num_modes <= solutions.size()));

    Hermes::vector<MeshFunctionSharedPtr<double> > slns;
    std::vector<std::vector<Function1D> > params;

    int num_parameters = parameters.size();
    for(int param_idx = 0; param_idx < num_parameters; param_idx++)
        params.push_back(std::vector<Function1D>());

    for(int mode_idx = 0; mode_idx < num_modes; mode_idx++){
        slns.push_back(solutions.at(mode_idx));

        for(int param_idx = 0; param_idx < num_parameters; param_idx++)
        {
            params[param_idx].push_back(parameters[param_idx][mode_idx]);
        }
    }

    slns.push_back(dirichlet_lift);

    CombinationFilter *cf = new CombinationFilter(slns, params, parameter_value);
    return MeshFunctionSharedPtr<double> (cf);
}


Function1D::~Function1D()
{
    if(points)
    {
        assert(values);
        delete[] values;
        delete[] points;
    }
    else
        assert(values == nullptr);
}
