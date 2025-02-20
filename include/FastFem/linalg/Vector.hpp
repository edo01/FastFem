#ifndef FASTFEM_VECTOR_HPP
#define FASTFEM_VECTOR_HPP

#include <vector>
#include <iostream>
#include <cmath>
#include <initializer_list>
#include <algorithm>

#include "FastFem/types/CommonTypes.hpp"

namespace fastfem{
namespace linalg{

using types::ff_index;

class Vector {
private:
    std::vector<double> data;

public:
    // Constructors
    Vector() = default;
    explicit Vector(ff_index size, double value = 0.0) : data(size, value) {}
    Vector(std::initializer_list<double> list) : data(list) {}
    
    // Utils
    inline ff_index size() const { return data.size(); }
    inline double& operator[](ff_index index) { return data[index]; }
    inline const double& operator[](ff_index index) const { return data[index]; }
    inline void fill(double value) { std::fill(data.begin(), data.end(), value); }
    inline double max() const { return *std::max_element(data.begin(), data.end()); }

    // Arithmetic Operations
    Vector operator+(const Vector& other) const;
    Vector operator-(const Vector& other) const;
    Vector operator*(double scalar) const;
    double dot(const Vector& other) const;
    double norm() const;
    static void axpby(double a, const Vector& x, double b, Vector& y);  // y = a*x + b*y
    static void axpy(double a, const Vector& x, Vector& y);  // y = a*x + y

    // I/O
    void print(std::string name = "") const;
};

} // namespace linalg
} // namespace FastFem

#endif // FASTFEM_VECTOR_HPP