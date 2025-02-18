#include "FastFem/mesh/MeshIO.hpp"

namespace fastfem{
namespace mesh{

template<unsigned int dim, unsigned int spacedim>
void MeshIO<dim,spacedim>::save_vtu(const std::string &filename) const{
    std::ofstream file;
    file.open(filename);
    file << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">\n";
    file << "<UnstructuredGrid>\n";
    file << "<Piece NumberOfPoints=\"" << mesh.vtx_count() << "\" NumberOfCells=\"" << mesh.elem_count() << "\">\n";
    file << "<Points>\n";
    file << "<DataArray type=\"Float64\" NumberOfComponents=\"" << 3 << "\" format=\"ascii\">\n"; 
    for (auto v = mesh.vtx_begin(); v != mesh.vtx_end(); ++v) {
        for(unsigned int i = 0; i < spacedim; i++)
            file << v->point.coords[i] << " ";
        if (spacedim < 3) {
            for(unsigned int i = spacedim; i < 3; i++)
                file << "0 ";
        }
        file << "\n";
    }
    file << "</DataArray>\n";
    file << "</Points>\n";
    file << "<Cells>\n";
    file << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (auto e = mesh.elem_begin(); e != mesh.elem_end(); ++e) {
        for (int i = 0; i < e->vertex_count(); ++i) {
            file << e->get_vertex(i) << " ";
        }
        file << "\n";
    }
    file << "</DataArray>\n";
    file << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    size_t offset = 0;
    for (auto e = mesh.elem_begin(); e != mesh.elem_end(); ++e) {
        offset += e->vertex_count();
        file << offset << "\n";
    }
    file << "</DataArray>\n";
    file << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (auto e = mesh.elem_begin(); e != mesh.elem_end(); ++e) {
        file << "5\n"; // 5 corresponds to VTK_TRIANGLE in VTK file format
    }
    file << "</DataArray>\n";
    file << "</Cells>\n";
    file << "</Piece>\n";
    file << "</UnstructuredGrid>\n";
    file << "</VTKFile>\n";
    file.close();
}

template<unsigned int dim, unsigned int spacedim>
void MeshIO<dim,spacedim>::save_msh(const std::string &filename) const{
    std::ofstream file;
    file.open(filename);
    file << "$MeshFormat\n";
    file << "2.2 0 8\n";
    file << "$EndMeshFormat\n";
    file << "$Nodes\n";
    file << mesh.vtx_count() << "\n";
    int id = 0;
    for (auto v = mesh.vtx_begin(); v != mesh.vtx_end(); ++v) {
        file << id++ << " ";
        for(unsigned int i = 0; i < spacedim; i++)
            file << v->point.coords[i] << " ";
        if (spacedim < 3) {
            for(unsigned int i = spacedim; i < 3; i++)
                file << "0 ";
        }
        file << "\n";
    }
    file << "$EndNodes\n";
    file << "$Elements\n";
    file << mesh.elem_count() << "\n";
    id=0;
    for (auto e = mesh.elem_begin(); e != mesh.elem_end(); ++e) {
        file << id++ << " ";
        file << "2 "; // 2 corresponds to VTK_TRIANGLE in Gmsh file format
        file << "3 "; // 3 corresponds to number of tags
        file << "1 1 1 "; // tags
        for (int i = 0; i < e->vertex_count(); ++i) {
            file << e->get_vertex(i) << " ";
        }
        file << "\n";
    }
    file << "$EndElements\n";
    file.close();
}


// Explicit instantiation
template class MeshIO<1,1>;
template class MeshIO<1,2>;
template class MeshIO<1,3>;
template class MeshIO<2,2>;
template class MeshIO<2,3>;
template class MeshIO<3,3>;


} // namespace mesh
} // namespace fastfem