#include "discretization.h"
#include "definitions.h"

MeshFunctionSharedPtr<double> PGDSolutions::get_filter(double parameter_value, int num_modes)
{
    if(num_modes == -1)
        num_modes = solutions.size();
    assert((num_modes >= 0) && (num_modes <= solutions.size()));

    Hermes::vector<MeshFunctionSharedPtr<double> > slns;
    std::vector<Function1D> params;
    for(int i = 0; i < num_modes; i++){
        slns.push_back(solutions.at(i));
        params.push_back(parameters.at(i));
    }

    slns.push_back(dirichlet_lift);

    CombinationFilter *cf = new CombinationFilter(slns, params, parameter_value);
    return MeshFunctionSharedPtr<double> (cf);
}


