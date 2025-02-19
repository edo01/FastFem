/**
 * @file DofHandler.hpp
 * 
 * The idea is the following: we need a class that is able to handle the degrees of freedom (DoFs) of a finite element space. 
 * This class should be able to distribute the DoFs on the mesh, to extract the DoFs on the boundary of the domain and to provide
 * an iterator that allows to loop over the cells of the mesh and to compute the local matrix and vector of the linear system.
 * 
 * Before going into the details of the implementation, let's see how the DoFHandler class is used in the MinSurFEM code:
 * We want avoid to have a double loop in the assembly of the linear system O(n_dofs^2), since most of these couples will be zero.
 * 
 * 
 * So the dof handler must expose a method that allows to iterate over the cells of the mesh, compute the local matrix and vector, 
 * and then add them to the global matrix and vector. Something like deal.II does:
 * for (auto &cell : dof_handler.cell_iterators())
 * {
 *  //  1. fill a full matrix of size dofs_per_cell x dofs_per_cell and a vector of size dofs_per_cell
 *  //    The difference with the deal.II code is that we do not provide a class that computes the values of the basis functions
 *  //    and their gradients at the quadrature points. The user should compute them and build the local matrix. 
 * 
 *  //  2. get the global indices of the DoFs of the current cell
 *  //  3. add the local matrix and vector to the global matrix and vector
 * }
 * 
 * In order to do this, the DoFHandler class must store for each cell the global indices of the DoFs of that cell.
 * To do this, we must remember that each DOF may be associated with a vertex, an edge, a face, or a cell of the mesh. A DOF on the vertex has support on the cells that contain that vertex, and so on. 
 * 
 * Where the information about the DoFs in each single cell is stored? Given a cell, we must be able to retrieve the points of the mesh that are associated with the DoFs of that cell. 
 *
 */



// numbering of degrees of freedom (DoFs):
// - we use an hash table with key the element and the value the array of the global indices of the DoFs of that element.
//      (element, [global_index_1, -1, -1, global_index_2, -1, -1, ...])
// - we iterate over all cells of the mesh:
//   - We get the indices of the DoFs of the current cell, and we set only the indices that are not already set
//   - As we set the indices, we increment a counter that is the number of DoFs of the mesh and we must update all the other elements of the hash table that share the same DoFs

// the last point is the most difficult to implement, since we must be able to retrieve all the elements that share the same DoFs.
// The idea is that, depending on the type of the element, the DoFs are placed on vertices, edges, faces, or the cell itself. For example, P1 elements have one DoF per vertex, P2 elements have also dof on the edges, P3 on the interior etc... . 
// One idea can be to store a VTAdjacency object that stores for each element its adjacent elements. When we set the indices of the DoFs of an element, we iterate over all the adjacent elements and we compute their intersection which is a simplex. We then compute the simplex that contains the updated DoFs and we check if the two simplices are the same. 

// We must add two methods to the simplex class: one that returns the barycentric coordinates of a point with respect to the vertices of the simplex.

#ifndef DOFHANDLER_HPP
#define DOFHANDLER_HPP

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <memory>
#include <algorithm>
#include <cassert>

#include "FastFem/mesh/Mesh.hpp"
#include "FastFem/types/CommonTypes.hpp"
#include "FastFem/fe/FESimplexP.hpp"

namespace fastfem{
namespace dof{

template <unsigned int dim, unsigned int spacedim = dim>
class DoFHandler
{
    static_assert(dim > 0, "The dimension must be greater than 0");
    static_assert(dim <= spacedim, "The dimension must be less or equal to the space it lives in");
    static_assert(spacedim <= 3, "The space dimension must be less or equal to 3");

    using global_dof_index_t = fastfem::types::global_dof_index_t;
    using local_dof_index_t  = fastfem::types::local_dof_index_t;

public:
    DoFHandler(const mesh::Mesh<dim, spacedim> &mesh, std::unique_ptr<fe::FESimplexP<dim, spacedim>> fe);

    /**
     * Distribute the degrees of freedom on the mesh. We offer a default implementation that works for 
     * generic finite element spaces.
     * 
     * The distribution of the DoFs is done in the following way: all the simplices of dimension 0 (vertices)
     * are visited. For each vertex, we assign a number of DoFs equal to the number of DoFs per vertex of
     * the finite element space. We then visit all the simplices of dimension 1 (edges) and we assign a number
     * of DoFs equal to the number of DoFs per edge of the finite element space. We then visit all the simplices
     * of dimension 2 (faces) and finally all the simplices of dimension 3 (tetrahedra). 
     * 
     * The Dofs are uniquely identified by their position inside a list local to the simplex they belong to. Each
     * simplex is identified by the couple in ascending order of the global indices of its vertices.
     * 
     */
    unsigned int distribute_dofs();

    /**
     * Get the global indices of the DoFs of the given element. The ordering of the DoFs is not coherent 
     * with the ordering of the local DoFs of the element. 
     */
    std::vector<global_dof_index_t> get_unordered_dofs_on_element(const mesh::MeshSimplex<dim, spacedim> &T) const;
    
    std::vector<global_dof_index_t> get_unordered_dofs_on_boundary(const mesh::MeshSimplex<dim-1, spacedim> &T) const;

    /**
     * Get the global indices of the DoFs of the given element. The ordering of the DoFs is coherent 
     * with the ordering of the local DoFs of the element. 
     */
    std::vector<global_dof_index_t> get_ordered_dofs_on_element(const mesh::MeshSimplex<dim, spacedim> &T) const;

    inline auto elem_begin() const { return mesh.elem_begin(); }
    inline auto elem_end() const { return mesh.elem_end(); }

    inline auto boundary_dofs_begin(size_t tag) const { return map_boundary_dofs.at(tag).begin(); }
    inline auto boundary_dofs_end(size_t tag) const { return map_boundary_dofs.at(tag).end(); }
    
    unsigned int get_n_dofs() const { return n_dofs; }

    inline size_t get_n_elements() const { return mesh.elem_count(); }

    void print_dofs() const;

private:
    const mesh::Mesh<dim, spacedim> &mesh;
    const std::unique_ptr<fe::FESimplexP<dim, spacedim>> fe;

    fastfem::types::global_dof_table<0> vertex_dofs;
    fastfem::types::global_dof_table<1> edge_dofs;
    fastfem::types::global_dof_table<2> face_dofs;
    fastfem::types::global_dof_table<3> cell_dofs;

    unsigned int n_dofs;

    std::unordered_map<size_t, std::unordered_set<global_dof_index_t>> map_boundary_dofs;

};

} // namespace dof
} // namespace fastfem

#endif // DOFHANDLER_HPP