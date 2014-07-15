#include "problem.h"
#include "definitions.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;

// This example shows how to solve a simple PDE that describes stationary
// heat transfer in an object consisting of two materials (aluminum and
// copper). The object is heated by constant volumetric heat sources
// (generated, for example, by a DC electric current). The temperature
// on the boundary is fixed. We will learn how to:
//
//   - load the mesh,
//   - perform initial refinements,
//   - create a H1 space over the mesh,
//   - define weak formulation,
//   - initialize matrix solver,
//   - assemble and solve the matrix system,
//   - output the solution and element orders in VTK format
//     (to be visualized, e.g., using Paraview),
//   - visualize the solution using Hermes' native OpenGL-based functionality.
//
// PDE: Poisson equation -div(LAMBDA grad u) - VOLUME_HEAT_SRC = 0.
//
// Boundary conditions: Dirichlet u(x, y) = FIXED_BDY_TEMP on the boundary.
//
// Geometry: L-Shape domain (see file domain.mesh).
//
// The following parameters can be changed:

const bool HERMES_VISUALIZATION = true;   // Set to "false" to suppress Hermes OpenGL visualization.
const bool VTK_VISUALIZATION = true;     // Set to "true" to enable VTK output.
const int P_INIT = 3;                     // Uniform polynomial degree of mesh elements.
const int INIT_REF_NUM = 1;               // Number of initial uniform mesh refinements.

// Problem parameters.




int main(int argc, char* argv[])
{
    // Load the mesh.
    MeshSharedPtr mesh(new Mesh);
    Hermes::Hermes2D::MeshReaderH2DXML mloader;
    mloader.load("mesh_cidlo.msh", mesh);

    // Refine all elements, do it INIT_REF_NUM-times.
    for (unsigned int i = 0; i < INIT_REF_NUM; i++)
        mesh->refine_all_elements();


    Hermes::vector<int> profile_coarse(0,0,0,0,0,0,0,0,0,5);
    int active_electrode = 7;
    double eps_rel_material = 5;

    ProblemConfiguration configuration(profile_coarse, active_electrode);
    ProblemDefinition_1 definition(eps_rel_material, configuration);

    // Initialize essential boundary conditions.
    Hermes::Hermes2D::DefaultEssentialBCConst<double> bc_essential_ground(definition.bc_ground, 0);
    Hermes::Hermes2D::DefaultEssentialBCConst<double> bc_essential_potential(definition.bc_potential, definition.POTENTIAL);
    Hermes::Hermes2D::EssentialBCs<double> bcs(Hermes::vector<EssentialBoundaryCondition<double> *> (&bc_essential_ground, &bc_essential_potential));

    // Initialize space.
    SpaceSharedPtr<double> space(new Hermes::Hermes2D::H1Space<double>(mesh, &bcs, P_INIT));
    SpaceSharedPtr<double> space_perm(new Hermes::Hermes2D::H1Space<double>(mesh, nullptr, 3));

    std::cout << "Ndofs: " << space->get_num_dofs() << std::endl;


    // Initialize the weak formulation.
    CustomWeakFormPoisson wf(definition);

    CustomWeakFormPermitivity wf_perm(definition);

    // Initialize the solution.
    MeshFunctionSharedPtr<double> sln(new Solution<double>);
    MeshFunctionSharedPtr<double> sln_perm(new Solution<double>);

    // Initialize linear solver.
    Hermes::Hermes2D::LinearSolver<double> linear_solver(&wf, space);
    Hermes::Hermes2D::LinearSolver<double> linear_solver_perm(&wf_perm, space_perm);

    // Solve the linear problem.
    try
    {
        linear_solver.solve();
        linear_solver_perm.solve();

        // Get the solution vector.
        double* sln_vector = linear_solver.get_sln_vector();
        double* sln_vector_perm = linear_solver_perm.get_sln_vector();

        // Translate the solution vector into the previously initialized Solution.
        Hermes::Hermes2D::Solution<double>::vector_to_solution(sln_vector, space, sln);
        Hermes::Hermes2D::Solution<double>::vector_to_solution(sln_vector_perm, space_perm, sln_perm);

        MyVolumetricIntegralCalculator integralCalculator(sln, definition);
        double* result = new double[1];
        result = integralCalculator.calculate(HERMES_ANY);
        std::cout << "Integral is " << result[0] << std::endl;

        // VTK output.
        if (VTK_VISUALIZATION)
        {
            // Output solution in VTK format.
            Hermes::Hermes2D::Views::Linearizer lin(FileExport);
            bool mode_3D = false;
            lin.save_solution_vtk(sln, "sln.vtk", "Temperature", mode_3D, 1);

            // Output mesh and element orders in VTK format.
            Hermes::Hermes2D::Views::Orderizer ord;
            ord.save_mesh_vtk(space, "mesh.vtk");
            ord.save_orders_vtk(space, "ord.vtk");
            ord.save_markers_vtk(space, "markers.vtk");
        }

        if (HERMES_VISUALIZATION)
        {
            // Visualize the solution.
            Hermes::Hermes2D::Views::ScalarView viewS("Solution", new Hermes::Hermes2D::Views::WinGeom(0, 0, 1500, 700));
            Hermes::Hermes2D::Views::ScalarView viewP("Permitivity", new Hermes::Hermes2D::Views::WinGeom(0, 700, 1500, 700));
            Hermes::Hermes2D::Views::OrderView viewSp("Space", new Hermes::Hermes2D::Views::WinGeom(0, 600, 1200, 600));
            viewSp.show(space);
            viewS.show(sln);
            viewP.show(sln_perm);
            viewS.wait_for_close();
        }
    }
    catch (Exceptions::Exception& e)
    {
        std::cout << e.info();
    }
    catch (std::exception& e)
    {
        std::cout << e.what();
    }
    return 0;
}
