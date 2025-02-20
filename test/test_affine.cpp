#include <iostream>
#include <cmath>
#include <memory>

#include "FastFem/mesh/MeshMaker.hpp"
#include "FastFem/fe/FESimplexP.hpp"
#include "FastFem/dof/DofHandler.hpp"

#include "FastFem/linalg/Vector.hpp"
#include "FastFem/mesh/MeshIO.hpp"
#include "FastFem/linalg/MatrixTools.hpp"


using namespace fastfem;

int main()
{
    auto f = [](double x, double y) { return  x*x -y*y; };

    mesh::SquareMaker mesh_maker(50);
    mesh::Mesh<2> mesh = mesh_maker.make_mesh();

    fe::FESimplexP2<2> fe;
    dof::DoFHandler<2> dof_handler(mesh);

    unsigned n_dofs = dof_handler.distribute_dofs(std::make_shared<fe::FESimplexP2<2>>(fe));

    linalg::Vector f_interpolated(n_dofs);

    linalg::MatrixTools::interpolate(f_interpolated, dof_handler, f);


    mesh::DataIO<2, 2> data_io(mesh, dof_handler, f_interpolated);
    data_io.save_vtx("interpolated.vtk");


        

    return 0;
}