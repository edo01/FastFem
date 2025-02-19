#ifndef MESHIO_HPP
#define MESHIO_HPP

#include <string>
#include <fstream>

#include "FastFem/mesh/Mesh.hpp"
#include "FastFem/dof/DofHandler.hpp"
#include "FastFem/linalg/Vector.hpp"

namespace fastfem {
namespace mesh {

using namespace fastfem::dof;
using namespace fastfem::linalg;

template <unsigned int dim, unsigned int spacedim = dim>
class MeshIO
{
public:
    MeshIO(Mesh<dim, spacedim> &mesh) : mesh(mesh) {}
    virtual ~MeshIO() = default;

    virtual void save_vtu(const std::string &filename) const;
    virtual void save_msh(const std::string &filename) const;

    //virtual Mesh<dim, spacedim> read(const std::string &filename) const = 0;

private:
    Mesh<dim, spacedim> &mesh;
};

/**
 * @brief Save the mesh in VTK format.
 * 
 * We use the VTK format to save the mesh. We plot the function
 */

template <unsigned int dim, unsigned int spacedim>
class DataIO
{
public:
    DataIO(Mesh<dim, spacedim> &mesh, DoFHandler<dim, spacedim> &dof_handler, Vector &solution);

    virtual ~DataIO() = default;

    virtual void save_vtx(const std::string &filename) const;

private:
    Mesh<dim, spacedim> &mesh;
    DoFHandler<dim, spacedim> &dof_handler;
    Vector &solution;
};

} // namespace mesh
} // namespace fastfem

#endif // MESHIO_HPP