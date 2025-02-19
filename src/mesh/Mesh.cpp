#include "FastFem/mesh/Mesh.hpp"

namespace fastfem{
namespace mesh{


using global_vertex_id = fastfem::types::global_vertex_id;
using global_edge_id = fastfem::types::global_edge_id;
using global_face_id = fastfem::types::global_face_id;
using global_cell_id = fastfem::types::global_cell_id;

/*
 * ----------------------VECTOR----------------------------
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
MeshSimplex<dim, spacedim>::MeshSimplex(const size_t v[dim + 1])     
{
    std::copy(v, v + dim + 1, vertices.begin());
}

template<unsigned int dim, unsigned int spacedim>
const std::array<global_vertex_id, MeshSimplex<dim, spacedim>::n_vertices> 
MeshSimplex<dim, spacedim>::get_vertex_indices() const {

    std::array<global_vertex_id, n_vertices> v_indices;
    unsigned int vertex_count = 0;

    for (size_t i = 0; i < dim + 1; ++i)
    {
        std::array<size_t, 1> v = {vertices[i]};
        v_indices[vertex_count++] = v;
    }
    return v_indices;
}

template<unsigned int dim, unsigned int spacedim>
const std::array<global_edge_id, MeshSimplex<dim, spacedim>::n_edges> 
MeshSimplex<dim, spacedim>::get_edges_indices() const {
    
    std::array<global_edge_id, n_edges> e_indices;
    unsigned int edge_count = 0;

    for (size_t i = 0; i < dim + 1; ++i)
    {
        for (size_t j = i + 1; j < dim + 1; ++j)
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
const std::array<global_face_id, MeshSimplex<dim, spacedim>::n_faces> 
MeshSimplex<dim, spacedim>::get_faces_indices() const {
    if(n_faces == 0) return {}; // in 1D there are no faces

    std::array<global_face_id, n_faces> f_indices;

    unsigned int face_count = 0;

    for (size_t i = 0; i < dim + 1; ++i){
        for (size_t j = i + 1; j < dim + 1; ++j){
            for (size_t k = j + 1; k < dim + 1; ++k){
                global_face_id face = {vertices[i], vertices[j], vertices[k]};
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
const std::array<global_cell_id, MeshSimplex<dim, spacedim>::n_cells> 
MeshSimplex<dim, spacedim>::get_cell_indices() const {
    if(n_cells == 0) return {}; // in 1D and 2D there are no cells
    
    std::array<global_cell_id, n_cells> c_indices;

    unsigned int cell_count = 0;

    for (size_t i = 0; i < dim + 1; ++i){
        for (size_t j = i + 1; j < dim + 1; ++j){
            for (size_t k = j + 1; k < dim + 1; ++k){
                for (size_t l = k + 1; l < dim + 1; ++l){
                   global_cell_id cell = {vertices[i], vertices[j], vertices[k], vertices[l]};
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
    for (unsigned int i = 0; i < dim + 1; ++i)
    {
        if (vertices[i] != s.vertices[i])
            return false;
    }
    return true;
}


/*
 * ----------------------MESH------------------------------
 */
// template<unsigned int dim, unsigned int spacedim>
// Simplex<dim, spacedim> Mesh<dim, spacedim>::get_Simplex(MeshSimplex<dim, spacedim> s) const
// {
//     Point<spacedim> p[dim + 1];
//     for (size_t i = 0; i < dim + 1; ++i)
//     {
//         p[i] = vertices[s.get_vertex(i)];
//     }
//     return Simplex<dim, spacedim>(p);
// }


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