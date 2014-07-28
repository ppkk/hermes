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

struct AllOnePerms : Perms
{
    AllOnePerms()
    {
        EPS_AIR = EPS0;
        EPS_EMPTY = EPS0;
        EPS_FULL = EPS0;
        EPS_KARTIT = EPS0;
    }
};


struct ProblemDefinition
{
    Hermes::vector<std::string> labels_air, labels_kartit, labels_full, labels_empty;
    Hermes::vector<std::string> bc_labels_ground, bc_labels_potential;
    double POTENTIAL;
    double SOURCE_TERM;
    char mesh_name[30];
    int P_INIT;                     // Uniform polynomial degree of mesh elements.
    int INIT_REF_NUM;               // Number of initial uniform mesh refinements.

};

struct ProblemDefinitionCidlo1 : ProblemDefinition
{
    ProblemDefinitionCidlo1(ProblemConfiguration configuration)
    {
        strcpy(mesh_name, "mesh_cidlo.msh");
        P_INIT = 2;
        INIT_REF_NUM = 0;

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

struct ProblemDefinitionUnitSquare : ProblemDefinition
{
    ProblemDefinitionUnitSquare(ProblemConfiguration configuration)
    {
        strcpy(mesh_name, "mesh_unit_square.msh");
        P_INIT = 4;
        INIT_REF_NUM = 3;

        POTENTIAL = 1.0;
        SOURCE_TERM = 0.0;

        // zemnici elektroda
        bc_labels_ground.push_back("Bottom");
        bc_labels_ground.push_back("Top");
        bc_labels_ground.push_back("Left");
        bc_labels_ground.push_back("Right");

        labels_full.push_back("Full");
    }
};

struct ProblemDefinitionUnitSquareDivided : ProblemDefinition
{
    ProblemDefinitionUnitSquareDivided(ProblemConfiguration configuration)
    {
        strcpy(mesh_name, "mesh_unit_square_divided.msh");
        P_INIT = 4;
        INIT_REF_NUM = 3;

        POTENTIAL = 1.0;
        SOURCE_TERM = 0.0;

        // zemnici elektroda
        bc_labels_ground.push_back("Bottom");
        bc_labels_ground.push_back("Top");
        bc_labels_ground.push_back("Left");
        bc_labels_ground.push_back("Right");

        labels_full.push_back("Full");
        labels_air.push_back("Air");
    }
};


#endif // PROBLEM_H
