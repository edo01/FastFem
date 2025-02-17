#ifndef FASTFEM_DOFHANDLER_HPP
#define FASTFEM_DOFHANDLER_HPP

#include <vector>
#include <random>


namespace fastfem{
namespace dof{

template <unsigned int dim, unsigned int spacedim = dim>
class DofHandler
{
private:
    std::vector<std::vector<unsigned int>> dofs_per_cell;

public:
    DofHandler(unsigned int n_elems)
    {
        dofs_per_cell.resize(n_elems);
        //rng
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<unsigned int> dis(0, 10);
        for(unsigned int i = 0; i < n_elems; ++i){
            unsigned int n_dofs = dis(gen);
            dofs_per_cell[i].resize(n_dofs);
            for(unsigned int j = 0; j < n_dofs; ++j){
                dofs_per_cell[i][j] = j;
            }
        }
    }
    
    inline unsigned int n_elements() const { return dofs_per_cell.size(); }
    inline const std::vector<unsigned int>& get_element_dofs(unsigned int i) const { return dofs_per_cell[i]; }
    
};

} // namespace dof
} // namespace fastfem

#endif // FASTFEM_DOFHANDLER_HPP