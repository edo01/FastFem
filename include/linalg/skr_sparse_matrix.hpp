#pragma once
#include<vector>
namespace linalg
{

struct SKRMatrix 
{
	size_t n_rows;
	size_t n_cols;
	size_t nnz;
    std::vector<double> data;
    std::vector<size_t> offset;
    std::vector<size_t> col_index;
}; 

}