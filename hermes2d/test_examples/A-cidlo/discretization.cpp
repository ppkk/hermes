#include "discretization.h"
#include "definitions.h"

MeshFunctionSharedPtr<double> PGDSolutions::get_filter(double parameter_value)
{
    Hermes::vector<MeshFunctionSharedPtr<double> > slns;
    for(int i = 0; i < solutions.size(); i++)
        slns.push_back(solutions.at(i));

    slns.push_back(dirichlet_lift);

    CombinationFilter *cf = new CombinationFilter(slns, parameters, parameter_value);
    return MeshFunctionSharedPtr<double> (cf);
}


