#ifndef FASTFEM_VECTOR_HPP
#define FASTFEM_VECTOR_HPP

#include <vector>
#include <iostream>
#include <cmath>
#include <initializer_list>
#include <algorithm>

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