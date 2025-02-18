#ifndef MESHIO_HPP
#define MESHIO_HPP

#include <string>
#include <fstream>

#include "FastFem/mesh/Mesh.hpp"

namespace fastfem {
namespace mesh {

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

} // namespace mesh
} // namespace fastfem

#endif // MESHIO_HPP