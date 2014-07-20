#ifndef PROBLEM_H
#define PROBLEM_H

#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;

const double EPS0 = 8.854e-12;

struct ProblemConfiguration
{
    ProblemConfiguration(Hermes::vector<int> profile_coarse, int active_electrode) : profile_coarse(profile_coarse), active_electrode(active_electrode) {}
    Hermes::vector<int> profile_coarse;
    int active_electrode;
};

struct Perms
{
    double EPS_AIR, EPS_KARTIT, EPS_FULL, EPS_EMPTY;

};

struct StandardPerms : Perms
{
    StandardPerms(double eps_rel_material)
    {
        EPS_AIR = EPS0;
        EPS_EMPTY = EPS0;
        EPS_FULL = eps_rel_material * EPS0;
        EPS_KARTIT = 4.5 * EPS0;
    }
};

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
        copy_from(origin);
    }

    Function1D& copy_from(const Function1D &origin)
    {
        n_intervals = origin.n_intervals;
        n_points = origin.n_points;
        bound_lo = origin.bound_lo;
        bound_hi = origin.bound_hi;
        points = new double[n_points];
        values = new double[n_points];
        for(int i = 0; i < n_points; i++)
        {
            points[i] = origin.points[i];
            values[i] = origin.values[i];
        }

        return *this;
    }


    ~Function1D()
    {
        if(points)
        {
            delete[] points;
            delete[] values;
        }
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
        for(int i = 0; i < n_intervals; i++)
        {
            assert(points[i] == ext.points[i]);
            result += (points[i]*values[i]*ext.values[i] + points[i+1]*values[i+1]*ext.values[i+1]) / 2.;
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

    double bound_lo, bound_hi;
    double *points, *values;
    int n_points, n_intervals;
};

struct ProblemDefinition
{
    Hermes::vector<std::string> labels_air, labels_kartit, labels_full, labels_empty;
    Hermes::vector<std::string> bc_labels_ground, bc_labels_potential;
    double POTENTIAL;
    double SOURCE_TERM;
    char mesh_name[30];
};

struct ProblemDefinition_1 : ProblemDefinition
{
    ProblemDefinition_1(ProblemConfiguration configuration)
    {
        strcpy(mesh_name, "mesh_cidlo.msh");

        POTENTIAL = 1.0;
        SOURCE_TERM = 0.0;

        int electrodes[8][4] = {{16, 17, 18, 19},
                                {20, 27, 28, 29},
                                {21, 30, 31, 32},
                                {22, 33, 34, 35},
                                {23, 36, 37, 38},
                                {24, 39, 40, 41},
                                {25, 42, 43, 44},
                                {26, 45, 46, 47}};



        // zemnici elektroda
        for(int i = 50; i <= 129; i++)
            bc_labels_ground.push_back(std::to_string(i));

        // klec
        for(int i = 10; i <= 15; i++)
            bc_labels_ground.push_back(std::to_string(i));

        for(int el_idx = 0; el_idx < 8; el_idx++)
        {
            for(int j = 0; j < 4; j++)
            {
                if(el_idx == configuration.active_electrode)
                    bc_labels_potential.push_back(std::to_string(electrodes[el_idx][j]));
                else
                    bc_labels_ground.push_back(std::to_string(electrodes[el_idx][j]));
            }
        }

        const int N_WIDTH = 80;
        const int N_HEIGHT = 18;
        int num_steps = configuration.profile_coarse.size();
        assert(N_WIDTH % num_steps == 0);
        int step_width = N_WIDTH / num_steps;

        Hermes::vector<int> profile;
        for(int i = 0; i < num_steps; i++)
            for(int j = 0; j < step_width; j++)
                profile.push_back(configuration.profile_coarse[i]);

        const int MAT_TO_LABEL_OFFSET = 11;
        for(int i = 0; i < N_WIDTH; i++)
        {
            for(int j = 0; j < N_HEIGHT; j++)
            {
                int idx = N_HEIGHT * i + j + MAT_TO_LABEL_OFFSET;
                if(profile[i] > j)
                    labels_full.push_back(std::to_string(idx));
                else
                {
                    //std::cout << idx << "\n";
                    labels_empty.push_back(std::to_string(idx));
                }
            }
        }

        labels_air.push_back("0");
        labels_air.push_back("10");
        labels_kartit.push_back("1");

    }
};

struct PGDSolutions
{
    PGDSolutions(ProblemDefinition definition, Perms perms, MeshSharedPtr mesh) : definition(definition), perms(perms), mesh(mesh) {}
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


    Function1D iteration_update_parameter();
    MeshFunctionSharedPtr<double> iteration_update_solution();

    void find_new_pair();

    std::vector<Function1D> parameters;
    std::vector<MeshFunctionSharedPtr<double> > solutions;

    Function1D actual_parameter;
    MeshFunctionSharedPtr<double> actual_solution;

    ProblemDefinition definition;
    Perms perms;
    MeshSharedPtr mesh;
};

#endif // PROBLEM_H
