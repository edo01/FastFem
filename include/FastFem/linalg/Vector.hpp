#ifndef VECTOR_HPP
#define VECTOR_HPP

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
    
    // Accessors
    std::size_t size() const { return data.size(); }
    double& operator[](std::size_t index) { return data[index]; }
    const double& operator[](std::size_t index) const { return data[index]; }

    // Arithmetic Operations
    Vector operator+(const Vector& other) const;
    Vector operator-(const Vector& other) const;
    Vector operator*(double scalar) const;
    Vector axpby(double a, const Vector& x, double b, const Vector& y) const;  // z = a*x + b*y
    void axpby(double a, const Vector& x, double b, Vector& y) const;   // y = a*x + b*y
    double dot(const Vector& other) const;
    double norm() const;

    // I/O
    void print() const;
};

} // namespace linalg
} // namespace FastFem

#endif // VECTOR_HPP