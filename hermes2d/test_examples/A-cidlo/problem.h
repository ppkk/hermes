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
    Function1D(double bound_lo, double bound_hi, int nintervals) : bound_lo(bound_lo), bound_hi(bound_hi)
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
    }

    ~Function1D()
    {
        delete[] points;
        delete[] values;
    }

    double int_length()
    {
        return (bound_hi - bound_lo) / n_intervals;
    }

    void print()
    {
        for(int i = 0; i < n_points; i++)
            std::cout << points[i] << "->" << values[i] << std::endl;
    }

    double value(double x)
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


    double int_F_F()
    {
        double result = 0;
        for(int i = 0; i < n_intervals; i++)
        {
            result += (values[i]*values[i] + values[i+1]*values[i+1]) / 2.;
        }

        return result * int_length();
    }

    double int_F_ExtF(Function1D ext)
    {
        double result = 0;
        for(int i = 0; i < n_intervals; i++)
        {
            assert(points[i] == ext.points[i]);
            result += (values[i]*ext.values[i] + values[i+1]*ext.values[i+1]) / 2.;
        }

        return result * int_length();
    }

    double int_x_F_F()
    {
        double result = 0;
        for(int i = 0; i < n_intervals; i++)
        {
            result += (points[i]*values[i]*values[i] + points[i+1]*values[i+1]*values[i+1]) / 2.;
        }

        return result * int_length();
    }

    double int_x_F_ExtF(Function1D ext)
    {
        double result = 0;
        for(int i = 0; i < n_intervals; i++)
        {
            assert(points[i] == ext.points[i]);
            result += (points[i]*values[i]*ext.values[i] + points[i+1]*values[i+1]*ext.values[i+1]) / 2.;
        }

        return result * int_length();
    }

    double bound_lo, bound_hi;
    double *points, *values;
    int n_points, n_intervals;
};

struct ProblemDefinition
{
    Hermes::vector<std::string> labels_air, labels_kartit, labels_full, labels_empty;
    Hermes::vector<std::string> bc_ground, bc_potential;
    double POTENTIAL;
    char mesh_name[30];
};

struct ProblemDefinition_1 : ProblemDefinition
{
    ProblemDefinition_1(ProblemConfiguration configuration)
    {
        strcpy(mesh_name, "mesh_cidlo.msh");

        POTENTIAL = 1.0;

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
            bc_ground.push_back(std::to_string(i));

        // klec
        for(int i = 10; i <= 15; i++)
            bc_ground.push_back(std::to_string(i));

        for(int el_idx = 0; el_idx < 8; el_idx++)
        {
            for(int j = 0; j < 4; j++)
            {
                if(el_idx == configuration.active_electrode)
                    bc_potential.push_back(std::to_string(electrodes[el_idx][j]));
                else
                    bc_ground.push_back(std::to_string(electrodes[el_idx][j]));
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

#endif // PROBLEM_H
