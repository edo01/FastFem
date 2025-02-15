#ifndef FASTFEM_VECTOR_HPP
#define FASTFEM_VECTOR_HPP

#include <vector>
#include <iostream>
#include <cmath>
#include <initializer_list>

namespace fastfem{
namespace linalg{

class Vector {
private:
    std::vector<double> data;

public:
    // Constructors
    Vector() = default;
    explicit Vector(std::size_t size, double value = 0.0) : data(size, value) {}
    Vector(std::initializer_list<double> list) : data(list) {}
    
    // Utils
    inline std::size_t size() const { return data.size(); }
    inline double& operator[](std::size_t index) { return data[index]; }
    inline const double& operator[](std::size_t index) const { return data[index]; }

    // Arithmetic Operations
    Vector operator+(const Vector& other) const;
    Vector operator-(const Vector& other) const;
    Vector operator*(double scalar) const;
    void axpby(double a, double b, const Vector& y);  // x = a*x + b*y
    void axpy(double a, const Vector& y);  // x = a*x + y
    double dot(const Vector& other) const;
    double norm() const;

    // I/O
    void print() const;
};

} // namespace linalg
} // namespace FastFem

#endif // FASTFEM_VECTOR_HPP