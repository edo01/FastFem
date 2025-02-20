#include "FastFem/mesh/Mesh.hpp"

namespace fastfem{
namespace mesh{

using global_vertex_index  = fastfem::types::global_vertex_index;
using global_edge_index    = fastfem::types::global_edge_index;
using global_face_index    = fastfem::types::global_face_index;
using global_cell_index    = fastfem::types::global_cell_index;
using local_vertex_index   = fastfem::types::local_vertex_index;

/*
 * ----------------------VERTEX----------------------------
 */
template<unsigned int spacedim>
bool Vertex<spacedim>::operator==(const Vertex<spacedim> &v) const 
{
    for (unsigned int i = 0; i < spacedim; ++i)
    {
        if (point.coords[i] != v.point.coords[i])
            return false;
    }
    return true;
}

/*
 * --------------------MESH SIMPLEX-------------------------
 */

template<unsigned int dim, unsigned int spacedim>
MeshSimplex<dim, spacedim>::MeshSimplex(const global_vertex_index v[dim + 1])     
{
    std::copy(v, v + dim + 1, vertices.begin());
}

template<unsigned int dim, unsigned int spacedim>
const std::array<global_vertex_index, MeshSimplex<dim, spacedim>::n_vertices> 
MeshSimplex<dim, spacedim>::get_vertex_indices() const {

    std::array<global_vertex_index, n_vertices> v_indices;
    std::copy(vertices.begin(), vertices.end(), v_indices.begin());

    return v_indices;
}

template<unsigned int dim, unsigned int spacedim>
const std::array<global_edge_index, MeshSimplex<dim, spacedim>::n_edges> 
MeshSimplex<dim, spacedim>::get_edges_indices() const {
    
    std::array<global_edge_index, n_edges> e_indices;
    unsigned int edge_count = 0;

    for (local_vertex_index i = 0; i < dim + 1; ++i)
    {
        for (local_vertex_index j = i + 1; j < dim + 1; ++j)
        {
            // always store the indices in ascending order
            if(vertices[i] < vertices[j])
                e_indices[edge_count++] = {vertices[i], vertices[j]};
            else
                e_indices[edge_count++] = {vertices[j], vertices[i]};
        }
    }
    
    assert(edge_count == n_edges);

    return e_indices;
}

template<unsigned int dim, unsigned int spacedim>
const std::array<global_face_index, MeshSimplex<dim, spacedim>::n_faces> 
MeshSimplex<dim, spacedim>::get_faces_indices() const {
    if(n_faces == 0) return {}; // in 1D there are no faces

    std::array<global_face_index, n_faces> f_indices;

    unsigned int face_count = 0;

    for (local_vertex_index i = 0; i < dim + 1; ++i){
        for (local_vertex_index j = i + 1; j < dim + 1; ++j){
            for (local_vertex_index k = j + 1; k < dim + 1; ++k){
                global_face_index face = {vertices[i], vertices[j], vertices[k]};
                std::sort(face.begin(), face.end());
                // copy starting from the face_count index
                f_indices[face_count++] = face;
            }
        }
    }

    assert(face_count == n_faces);

    return f_indices;
}

template<unsigned int dim, unsigned int spacedim>
const std::array<global_cell_index, MeshSimplex<dim, spacedim>::n_cells> 
MeshSimplex<dim, spacedim>::get_cell_indices() const {
    if(n_cells == 0) return {}; // in 1D and 2D there are no cells
    
    std::array<global_cell_index, n_cells> c_indices;

    unsigned int cell_count = 0;

    for (local_vertex_index i = 0; i < dim + 1; ++i){
        for (local_vertex_index j = i + 1; j < dim + 1; ++j){
            for (local_vertex_index k = j + 1; k < dim + 1; ++k){
                for (local_vertex_index l = k + 1; l < dim + 1; ++l){
                   global_cell_index cell = {vertices[i], vertices[j], vertices[k], vertices[l]};
                    std::sort(cell.begin(), cell.end());
                    c_indices[cell_count++] = cell;                        
                }
            }
        }
    }

    assert(cell_count == n_cells);

    return c_indices;
}

template<unsigned int dim, unsigned int spacedim>
bool MeshSimplex<dim, spacedim>::operator==(const MeshSimplex<dim, spacedim> &s) const
{
    for (local_vertex_index i = 0; i < dim + 1; ++i)
    {
        if (vertices[i] != s.vertices[i])
            return false;
    }
    return true;
}

/**
 * ---------------Template instantiations-------------------
 */

// vertex
template struct Vertex<1>;
template struct Vertex<2>;
template struct Vertex<3>;

// mesh simplex
template class MeshSimplex<0, 1>;
template class MeshSimplex<0, 2>;
template class MeshSimplex<0, 3>;
template class MeshSimplex<1, 1>;
template class MeshSimplex<1, 2>;
template class MeshSimplex<1, 3>;
template class MeshSimplex<2, 2>;
template class MeshSimplex<2, 3>;
template class MeshSimplex<3, 3>;

} // namespace mesh
} // namespace fastfem