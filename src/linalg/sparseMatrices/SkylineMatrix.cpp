#include "FastFem/linalg/sparseMatrices/SkylineMatrix.hpp"

namespace fastfem{
namespace linalg{

SkylineMatrix::SkylineMatrix(size_t n_cols, const SkylinePattern& skyline) :
  SparseMatrix(skyline.skyline_rows.size() - 1, n_cols),
    values(skyline.skyline_rows.back()),
    base_skyline(std::make_shared<SkylinePattern>(skyline))
    {
        if(n_cols == 0 || skyline.skyline_rows.size() < 2){
            throw std::invalid_argument("SkylineMatrix::SkylineMatrix(): invalid dimensions");
        }
        if(skyline.skyline_rows.size() - 1 != n_cols){
            throw std::invalid_argument("SkylineMatrix::SkylineMatrix(): matrix must be square");
        }
    }

const double &SkylineMatrix::get_entry(size_t i, size_t j) const
{
    if (j > i) {
        std::swap(i, j);
    }

    size_t row_start = (base_skyline->skyline_rows)[i];
    size_t row_end = (base_skyline->skyline_rows)[i + 1];
    size_t row_length = row_end - row_start;

    // Compute the first column stored in row i
    size_t first_col_in_row = i - row_length + 1;

    // Check if A(i, j) is an implicit zero
    if (j < first_col_in_row) {
        static double zero = 0.0;
        return zero;
    }

    // Find position of A(i, j) in values array
    size_t position = row_start + (j - first_col_in_row);
    return values[position]; 
}

Vector SkylineMatrix::gemv(const Vector& x) const
{
    size_t row_start, row_end, row_length, first_col_in_row, col_index;

    if (n_cols != x.size()) {
        throw std::invalid_argument("SkylineMatrix::gemv(): incompatible dimensions");
    }

    Vector y(n_rows);

    for (size_t i = 0; i < n_rows; ++i)
    {
        //row_start = skyline->at(i);
        //row_end = skyline->at(i + 1);
        row_start = base_skyline->skyline_rows[i];
        row_end = base_skyline->skyline_rows[i + 1];
        row_length = row_end - row_start;

        // First column stored in row i
        first_col_in_row = i - row_length + 1;

        col_index = first_col_in_row;
        for (size_t k = row_start; k < row_end; ++k)
        {
            y[i] += values[k] * x[col_index];

            // Since the matrix is symmetric, update y[col_index] (except for diagonal elements)
            if (i != col_index) {
                y[col_index] += values[k] * x[i];
            }

            ++col_index;
        }
    }

    return y;
}

//TO BE CANCELLED
void SkylineMatrix::set_values(const std::vector<double>& values)
{
    if (values.size() != this->values.size()) {
        throw std::invalid_argument("SkylineMatrix::set_values(): incompatible dimensions");
    }
    this->values = values;
}

void SkylineMatrix::set_row_col_to_zero(size_t i)
{
    size_t row_start = (base_skyline->skyline_rows)[i];
    size_t row_end = (base_skyline->skyline_rows)[i + 1];

    #pragma omp parallel
    {
        #pragma omp for
        for (size_t k = row_start; k < row_end; ++k)
        {
            values[k] = 0.0;
        }

        #pragma omp for
        for(size_t j = 0; j < n_rows; ++j){
            set_entry_to_zero(i, j);
        }
    }

}

void SkylineMatrix::print_pattern() const
{
    for (size_t i = 0; i < n_rows; ++i)
    {
        for (size_t j = 0; j <= i; ++j)
        {
            bool found = false;
            for (size_t k = (base_skyline->skyline_rows)[i]; k < (base_skyline->skyline_rows)[i + 1]; ++k)
            {
                if (j == i - ((base_skyline->skyline_rows)[i + 1] - k - 1))  
                {
                    std::cout << "1 ";
                    found = true;
                    break;
                }
            }
            if (!found)
            {
                std::cout << "0 ";
            }
        }
        std::cout << std::endl;

    }

    std::cout << "\nSkyline rows: ";
    for(size_t i=0; i<this->base_skyline->skyline_rows.size(); i++){
        std::cout << this->base_skyline->skyline_rows[i] << " ";
    }
    std::cout << std::endl;
    std::cout << "Values: ";
    for(size_t i=0; i<this->values.size(); i++){
        std::cout << this->values[i] << " ";
    }
    std::cout << std::endl;
}

/**
 * @brief Computes the in-place Cholesky factorization: A = LL^T
 */
void SkylineMatrix::cholesky_factorize() {
    for (size_t i = 0; i < n_rows; ++i) {
        size_t row_start = (base_skyline->skyline_rows)[i];
        size_t row_end = (base_skyline->skyline_rows)[i + 1];
        size_t first_col = i - (row_end - row_start) + 1;

        // Compute L(i, i)
        double sum = 0.0;
        for (size_t k = first_col; k < i; ++k) {
            sum += (*this)(i, k) * (*this)(i, k);
        }

        values[row_end - 1] = std::sqrt((*this)(i, i) - sum);

        for (size_t j = i + 1; j < n_rows; ++j) {
            size_t row_start_j = (base_skyline->skyline_rows)[j];
            size_t row_end_j = (base_skyline->skyline_rows)[j + 1];
            size_t row_length_j = row_end_j - row_start_j;  // Length of stored values in row j
            size_t first_col_j = j - row_length_j + 1;
        
            // Stop looping if column i is outside the skyline storage of row j
            if (first_col_j > i) continue;
        
            double sum = 0.0;
            for (size_t k = first_col; k < i; ++k) {
                sum += (*this)(j, k) * (*this)(i, k);
            }
        
            values[row_start_j + (i - first_col_j)] = ((*this)(j, i) - sum) / values[row_end - 1];
        }        
    }
}


/**
 * @brief Solves Ax = b using the precomputed Cholesky factor L.
 */
Vector SkylineMatrix::cholesky_solve(const Vector& b) const {
    if (n_rows != b.size()) {
        throw std::invalid_argument("Incompatible dimensions");
    }

    Vector x = b;

    // Forward Substitution: Solve L * y = b
    for (size_t i = 0; i < n_rows; ++i) {
        size_t row_start = (base_skyline->skyline_rows)[i];
        size_t row_end = (base_skyline->skyline_rows)[i + 1];
        size_t first_col = i - (row_end - row_start) + 1;

        for (size_t k = row_start; k < row_end - 1; ++k) {
            x[i] -= values[k] * x[first_col + (k - row_start)];
        }
        x[i] /= values[row_end - 1];
    }

    // Backward Substitution: Solve L^T * x = y
    for (size_t i = n_rows; i > 0; --i) {       
        size_t row_start = (base_skyline->skyline_rows)[i - 1];
        size_t row_end = (base_skyline->skyline_rows)[i];
        size_t first_col = i - (row_end - row_start);

        x[i - 1] /= values[row_end - 1];

        for (size_t k = row_start; k < row_end - 1; ++k) {
            x[first_col + (k - row_start)] -= values[k] * x[i - 1];
        }
    }

    return x;
}

void SkylineMatrix::set_entry(size_t i, size_t j, double value) {
    if (i < j) {
        std::swap(i, j);
    }

    size_t row_start = (base_skyline->skyline_rows)[i];
    size_t row_end = (base_skyline->skyline_rows)[i + 1];
    size_t row_length = row_end - row_start;
    size_t first_col = i - row_length + 1;

    if (j < first_col) {
        throw std::out_of_range("SkylineMatrix::insert_entry(): Position is outside skyline storage.");
    }

    // Compute the index in the values array
    size_t index = row_start + (j - first_col);

    values[index] = value;
}

void SkylineMatrix::accumulate_entry(size_t i, size_t j, double value) {
    if (i < j) {
        std::swap(i, j);
    }

    size_t row_start = (base_skyline->skyline_rows)[i];
    size_t row_end = (base_skyline->skyline_rows)[i + 1];
    size_t row_length = row_end - row_start;
    size_t first_col = i - row_length + 1;

    // Compute the index in the values array
    size_t index = row_start + (j - first_col);

    values[index] += value;
}

void SkylineMatrix::set_entry_to_zero(size_t i, size_t j)
{
    if (i < j) {
        std::swap(i, j);
    }

    size_t row_start = (base_skyline->skyline_rows)[i];
    size_t row_end = (base_skyline->skyline_rows)[i + 1];
    size_t row_length = row_end - row_start;
    size_t first_col = i - row_length + 1;

    if (j < first_col) {
        return;
    }

    // Compute the index in the values array
    size_t index = row_start + (j - first_col);

    values[index] = 0.0;
}


template <unsigned int dim, unsigned int spacedim>
SkylinePattern SkylinePattern::create_from_dof_handler(const fastfem::dof::DoFHandler<dim, spacedim>& dof_handler)
{
    std::vector<std::set<unsigned int>> dof_interactions(dof_handler.get_n_dofs());
    
    for (auto it = dof_handler.elem_begin(); it != dof_handler.elem_end(); ++it) {
        const auto& elem = *it;

        const auto& dofs = dof_handler.get_ordered_dofs_on_element(elem);
        for(int j = 0; j < dofs.size(); ++j){
            for(int k = j; k < dofs.size(); ++k){
                dofs[j] > dofs[k] ? dof_interactions[dofs[j]].insert(dofs[k]) : dof_interactions[dofs[k]].insert(dofs[j]); 
            }
        }
    }

    std::vector<size_t> skyline_rows(dof_handler.get_n_dofs() + 1);

    for(unsigned int i = 0; i < dof_interactions.size(); ++i){
        skyline_rows[i + 1] = skyline_rows[i] + (i - *dof_interactions[i].begin() + 1);
    }

    return SkylinePattern(skyline_rows);
}

// instantiate the templates
template SkylinePattern SkylinePattern::create_from_dof_handler<1,1>(const fastfem::dof::DoFHandler<1,1>& dof_handler);
template SkylinePattern SkylinePattern::create_from_dof_handler<1,2>(const fastfem::dof::DoFHandler<1,2>& dof_handler);
template SkylinePattern SkylinePattern::create_from_dof_handler<1,3>(const fastfem::dof::DoFHandler<1,3>& dof_handler);
template SkylinePattern SkylinePattern::create_from_dof_handler<2,2>(const fastfem::dof::DoFHandler<2,2>& dof_handler);
template SkylinePattern SkylinePattern::create_from_dof_handler<2,3>(const fastfem::dof::DoFHandler<2,3>& dof_handler);
template SkylinePattern SkylinePattern::create_from_dof_handler<3,3>(const fastfem::dof::DoFHandler<3,3>& dof_handler);

} // namespace linalg
} // namespace FastFem