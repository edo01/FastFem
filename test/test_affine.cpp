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

        mesh::Simplex<2, 2> triangle = mesh.get_Simplex(elem);

        std::vector<types::global_dof_index> global_dofs = dof_handler.get_ordered_dofs_on_element(elem);

        for(types::local_dof_index i = 0; i < global_dofs.size(); ++i)
        {
            mesh::Point<2> dof_coords = fe.get_dof_coords(triangle, i);
            f_interpolated[global_dofs[i]] = f(dof_coords.coords[0], dof_coords.coords[1]);
        }
    }


    mesh::DataIO<2, 2> data_io(mesh, dof_handler, f_interpolated);
    data_io.save_vtx("interpolated.vtk");


        

    return 0;
}