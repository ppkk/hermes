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

struct ProblemDefinition
{
    Hermes::vector<std::string> labels_air, labels_kartit, labels_full, labels_empty;
    Hermes::vector<std::string> bc_ground, bc_potential;
    double EPS_AIR, EPS_KARTIT, EPS_FULL, EPS_EMPTY;
    double POTENTIAL;
    int ACTIVE_ELECTRODE;
};

struct ProblemDefinition_1 : ProblemDefinition
{
    ProblemDefinition_1(double eps_rel_material, ProblemConfiguration configuration)
    {
        EPS_AIR = EPS0;
        EPS_EMPTY = EPS0;
        EPS_FULL = eps_rel_material * EPS0;
        EPS_KARTIT = 4.5 * EPS0;

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


        // todo: refact
        Hermes::vector<int> profile;
        for(int i = 0; i < 10; i++)
            for(int j = 0; j < 8; j++)
                profile.push_back(configuration.profile_coarse[i]);

        const int N_WIDTH = 80;
        const int N_HEIGHT = 18;
        assert(N_WIDTH % configuration.profile_coarse.size() == 0);

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
