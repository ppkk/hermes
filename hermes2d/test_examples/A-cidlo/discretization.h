#ifndef DISCRETIZATION_H
#define DISCRETIZATION_H

#include "hermes2d.h"
#include "problem.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;

class CombinationFilter;

struct Function1D
{
    Function1D()
    {
        n_intervals = n_points = bound_hi = bound_lo = -12345;
        points = values = nullptr;
    }

    Function1D(double bound_lo, double bound_hi, int nintervals, double initial_value = 0) : bound_lo(bound_lo), bound_hi(bound_hi)
    {
        assert(bound_hi > bound_lo);
        n_intervals = nintervals;
        n_points = nintervals + 1;
        points = new double[n_points];
        values = new double[n_points];
        double int_len = int_length();
        points[0] = bound_lo;
        for(int i = 1; i < n_intervals ; i++)
            points[i] = points[i-1] + int_len;
        points[n_points - 1] = bound_hi;

        for(int i = 0; i < n_points; i++)
            values[i] = initial_value;
    }

    Function1D(const Function1D &origin)
    {
        copy_from(origin);
    }

    Function1D& operator=(const Function1D &origin)
    {
        return copy_from(origin);
    }

    Function1D& copy_from(const Function1D &origin)
    {
        n_intervals = origin.n_intervals;
        n_points = origin.n_points;
        bound_lo = origin.bound_lo;
        bound_hi = origin.bound_hi;
        if(n_points > 0)
        {
            points = new double[n_points];
            values = new double[n_points];
        }
        else
        {
            points = nullptr;
            values = nullptr;
        }

        for(int i = 0; i < n_points; i++)
        {
            points[i] = origin.points[i];
            values[i] = origin.values[i];
        }

        return *this;
    }


    ~Function1D();

    void normalize_first_to_one()
    {
        assert(values[0] != 0);
        double value_0 = values[0];
        for(int i = 0; i < n_points; i++)
            values[i] /= value_0;
    }

    double int_length() const
    {
        return (bound_hi - bound_lo) / n_intervals;
    }

    void print() const
    {
        for(int i = 0; i < n_points; i++)
            std::cout << points[i] << "->" << values[i] << std::endl;
    }

    double value(double x) const
    {
        assert(x >= bound_lo);
        assert(x <= bound_hi);

        // todo: improve
        int idx = 0;
        while(x > points[idx+1])
            idx++;
        assert(idx < n_intervals);
        return values[idx] + (values[idx+1] - values[idx]) * ( x -points[idx]) / (points[idx+1] - points[idx]);
    }

    void print_short() const
    {
        for(int i = 0; i < n_points; i++)
        {
            std::cout << values[i] << ", ";
        }
        std::cout << std::endl;
    }

    double int_F() const
    {
        double result = 0;
        for(int i = 0; i < n_intervals; i++)
        {
            result += (values[i] + values[i+1]) / 2.;
        }

        return result * int_length();
    }

    double int_x_F() const
    {
        double result = 0;
        for(int i = 0; i < n_intervals; i++)
        {
            result += (points[i] * values[i] + points[i+1] * values[i+1]) / 2.;
        }

        return result * int_length();
    }

    double int_F_F() const
    {
        double result = 0;
        for(int i = 0; i < n_intervals; i++)
        {
            result += (values[i]*values[i] + values[i+1]*values[i+1]) / 2.;
        }

        return result * int_length();
    }

    double int_F_ExtF(Function1D ext) const
    {
        double result = 0;
        for(int i = 0; i < n_intervals; i++)
        {
            assert(points[i] == ext.points[i]);
            result += (values[i]*ext.values[i] + values[i+1]*ext.values[i+1]) / 2.;
        }

        return result * int_length();
    }

    double int_x_F_F() const
    {
        double result = 0;
        for(int i = 0; i < n_intervals; i++)
        {
            result += (points[i]*values[i]*values[i] + points[i+1]*values[i+1]*values[i+1]) / 2.;
        }

        return result * int_length();
    }

    double int_x_F_ExtF(Function1D ext) const
    {
        double result = 0;
        assert(n_intervals == ext.n_intervals);
        for(int i = 0; i < n_intervals; i++)
        {
            assert(points[i] == ext.points[i]);
            result += (points[i]*values[i]*ext.values[i] + points[i+1]*values[i+1]*ext.values[i+1]) / 2.;
        }

        return result * int_length();
    }

    // dependence of permitivity on the value of parameter for given element
    // if parameter is below element given by element_height,  value is empty
    // if it is above, value is full
    // approximated linearly inside the element given by element_height (if noninteger values of parameter are allowed)
    static double perm_on_parameter(double parameter, int element_height, Perms perms)
    {
        assert((element_height >= 0) && (element_height < N_HEIGHT_COARSE));
        if(parameter < element_height)
            return perms.EPS_EMPTY;
        else if(parameter > element_height + 1)
            return perms.EPS_FULL;
        else
        {
            double ratio = parameter - element_height;
            assert((ratio >= 0) && (ratio <= 1));
            double result = perms.EPS_EMPTY + ratio * (perms.EPS_FULL - perms.EPS_EMPTY);
            assert((result >= perms.EPS_EMPTY) && (result <= perms.EPS_FULL));
            return result;
        }
    }

    // integrals for the case of columns with the peremitivity FULL
    // calculated for element_heights 0 ...N_HEIGHT_COARSE - 1
    // parameter has to be defined in [0, N_HEIGHT_COARSE]
    double int_epsx_F(int element_height, Perms perms) const
    {
        double result = 0;
        for(int i = 0; i < n_intervals; i++)
        {
            double perm_i = perm_on_parameter(points[i], element_height, perms);
            double perm_i_plus_1 = perm_on_parameter(points[i+1], element_height, perms);

            result += (perm_i * values[i] + perm_i_plus_1 * values[i+1]) / 2.;
        }

        return result * int_length();
    }

    double int_epsx_F_F(int element, Perms perms) const
    {
        double result = 0;
        for(int i = 0; i < n_intervals; i++)
        {
            double perm_i = perm_on_parameter(points[i], element, perms);
            double perm_i_plus_1 = perm_on_parameter(points[i+1], element, perms);

            result += (perm_i * values[i] * values[i] + perm_i_plus_1 * values[i+1] * values[i+1]) / 2.;
        }

        return result * int_length();
    }


    double diference(Function1D second) const
    {
        assert(second.n_points == this->n_points);

        double result = 0;
        for(int i = 0; i < n_points; i++)
        {
            result += sqr(values[i] - second.values[i]);
        }

        return Hermes::sqrt(result);
    }

    double norm() const
    {

        double result = 0;
        for(int i = 0; i < n_points; i++)
        {
            result += sqr(values[i]);
        }

        return Hermes::sqrt(result);
    }

    static void test()
    {
        Function1D fn(0., 1., 2);
        fn.values[0] = 1;
        fn.values[1] = 2;
        fn.values[2] = 1;
        fn.print();
        std::cout << fn.int_F_F() << std::endl;
        std::cout << fn.value(0.9) << std::endl;

        Function1D a(1.4,1.5,5);
        Function1D b;
        b = a;
        Function1D c(a);

        a.print();
        b.print();
        c.print();

        assert(a.int_F_ExtF(a) == a.int_F_F());
        assert(a.int_x_F_ExtF(a) == a.int_x_F_F());

        a.values[0] = 4;
        a.values[1] = 3;
        a.values[2] = 1;
        b.values[0] = -5;
        b.values[1] = 2;
        b.values[2] = 11;
        assert(a.int_x_F_ExtF(b) == b.int_x_F_ExtF(a));
        assert(a.int_F_ExtF(b) == b.int_F_ExtF(a));

        int N = 100;
        Function1D x(2., 3., N);
        for(int i = 0; i < x.n_points; i++)
        {
            x.values[i] = std::pow(x.points[i], 2);
        }
        int p = 5;
        double val_ref = (std::pow(3,p) - std::pow(2, p)) / p;
        std::cout << "should be similar " << x.int_F_F() << ", " << val_ref << std::endl;
        p = 6;
        val_ref = (std::pow(3,p) - std::pow(2, p)) / p;
        std::cout << "should be similar " << x.int_x_F_F() << ", " << val_ref << std::endl;

        std::cout << std::endl << "Testing columns integrals " << std::endl << std::endl;
        Perms perms;
        perms.EPS_EMPTY = 1;
        perms.EPS_FULL = 5;
        for(int element = 0; element < N_HEIGHT_COARSE; element++)
        {
            std::cout << "element " << element << ": ";
            for(double parameter = 0; parameter <= N_HEIGHT_COARSE; parameter += 0.5)
            {
                std::cout << perm_on_parameter(parameter, element, perms) << ", ";
            }
            std::cout << std::endl;
        }
    }

    double bound_lo, bound_hi;
    double *points, *values;
    int n_points, n_intervals;
};


struct PGDSolutions
{
    PGDSolutions(ProblemDefinition *definition, Perms perms, MeshSharedPtr mesh, int num_parameters) :
        definition(definition), perms(perms), mesh(mesh), num_parameters(num_parameters)
    {
        for(int i = 0; i < num_parameters; i++)
        {
            parameters.push_back(std::vector<Function1D>());
            Function1D fn;
            actual_parameter.push_back(fn);
        }
    }

    double get_pt_value(double x, double y, double eps, int use_modes = -1)
    {
        for(int i = 0; i < num_parameters; i++)
            assert(solutions.size() == parameters[i].size());

        if(use_modes == -1)
            use_modes = solutions.size();

        double result = 0;
        for(int i = 0; i < use_modes; i++)
        {
            double value = solutions.at(i)->get_pt_value(x, y)->val[0];
            for(int j = 0; j < num_parameters; j++)
            {
                value *= parameters[j].at(i).value(eps);
            }
            result += value;
        }

        if(definition->use_dirichlet_lift())
            result += dirichlet_lift->get_pt_value(x, y)->val[0];

        return result;
    }

    // if num_modes == -1 use all of them
    MeshFunctionSharedPtr<double> get_filter(double parameter_value, int num_modes = -1);

    Function1D iteration_update_parameter_changing_perm();
    MeshFunctionSharedPtr<double> iteration_update_solution_changing_perm();
    void find_new_pair();

    Function1D iteration_update_parameter_columns();

    int num_parameters;

    // array of arrays of modes, that correspond to individual parameters
    // parameters[parameter_idx][mode_idx]
    std::vector<std::vector<Function1D> > parameters;
    std::vector<MeshFunctionSharedPtr<double> > solutions;

    MeshFunctionSharedPtr<double> dirichlet_lift;

    // array of individual parameters
    // actual_parameter[parameter_idx]
    std::vector<Function1D> actual_parameter;
    MeshFunctionSharedPtr<double> actual_solution;

    ProblemDefinition *definition;
    Perms perms;
    MeshSharedPtr mesh;
};

#endif // DISCRETIZATION_H
