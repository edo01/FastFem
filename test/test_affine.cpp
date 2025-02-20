#include "FastFem/mesh/MeshMaker.hpp"
#include "FastFem/fe/FESimplexP.hpp"
#include "FastFem/dof/DofHandler.hpp"

#include "FastFem/linalg/Vector.hpp"
#include "FastFem/mesh/MeshIO.hpp"

using namespace fastfem;

int main()
{
    auto f = [](double x, double y) { return std::sqrt(1 - x*x -y*y); };

    mesh::SquareMaker mesh_maker(50);
    mesh::Mesh<2> mesh = mesh_maker.make_mesh();

    fe::FESimplexP2<2> fe;
    dof::DoFHandler<2> dof_handler(mesh, std::make_unique<fe::FESimplexP2<2>>());

    unsigned n_dofs = dof_handler.distribute_dofs();

    linalg::Vector f_interpolated(n_dofs);

    for(auto it = dof_handler.elem_begin(); it != dof_handler.elem_end(); ++it)
    {
        auto &elem = *it;
        std::cout << "Element: (" << elem.get_vertex(0) << ", " << elem.get_vertex(1) << ", " << elem.get_vertex(2) << ")" << std::endl;

        mesh::Simplex<2, 2> triangle = mesh.get_Simplex(elem);
        std::cout << "Vertices: (" << triangle.get_vertex(0).coords[0] << ", " << triangle.get_vertex(0).coords[1] << "), ";
        std::cout << "(" << triangle.get_vertex(1).coords[0] << ", " << triangle.get_vertex(1).coords[1] << "), ";
        std::cout << "(" << triangle.get_vertex(2).coords[0] << ", " << triangle.get_vertex(2).coords[1] << ")" << std::endl;

        std::vector<types::global_dof_index_t> global_dofs = dof_handler.get_ordered_dofs_on_element(elem);
        std::cout << "Global ordered dofs: ";
        for(auto dof : global_dofs)
        {
            std::cout << dof << ", ";
        }
        std::cout << std::endl;
        std::cout << "dof_coords: ";
        for(types::local_dof_index_t i = 0; i < global_dofs.size(); ++i)
        {
            mesh::Point<2> dof_coords = fe.get_dof_coords(triangle, i);
            std::cout << "(" << dof_coords.coords[0] << ", " << dof_coords.coords[1] << "), ";
            f_interpolated[global_dofs[i]] = f(dof_coords.coords[0], dof_coords.coords[1]);
        }
        std::cout << std::endl << std::endl;
    }

    f_interpolated.print("f_interpolated");
    std::cout << "Max of interpolated function: " << f_interpolated.max() << std::endl;

    mesh::DataIO<2, 2> data_io(mesh, dof_handler, f_interpolated);
    data_io.save_vtx("interpolated.vtk");


        

    return 0;
}