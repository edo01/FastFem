#include "FastFem/linalg/Vector.hpp"

namespace fastfem {
namespace linalg {


Vector Vector::operator+(const Vector& other) const {
    Vector result(data.size());
    #pragma omp parallel for
    for (std::size_t i = 0; i < data.size(); ++i) {
        result.data[i] = data[i] + other.data[i];
    }
    return result;
}

Vector Vector::operator-(const Vector& other) const {
    Vector result(data.size());
    #pragma omp parallel for
    for (std::size_t i = 0; i < data.size(); ++i) {
        result.data[i] = data[i] - other.data[i];
    }
    return result;
}

Vector Vector::operator*(double scalar) const {
    Vector result(data.size());
    #pragma omp parallel for
    for (std::size_t i = 0; i < data.size(); ++i) {
        result.data[i] = data[i] * scalar;
    }
    return result;
}

// y = a*x + b*y
void Vector::axpby(double a, const Vector& x, double b, Vector& y)
{
    #pragma omp parallel for
    for (std::size_t i = 0; i < x.size(); ++i) {
        y[i] = a * x[i] + b * y[i];
    }
}

// y = a*x + y
void Vector::axpy(double a, const Vector& x, Vector& y)
{
    #pragma omp parallel for
    for (std::size_t i = 0; i < x.size(); ++i) {
        y[i] = a * x[i] + y[i];
    }
}

double Vector::dot(const Vector& other) const {
    double result = 0.0;
    #pragma omp parallel for reduction(+:result)
    for (std::size_t i = 0; i < data.size(); ++i) {
        result += data[i] * other.data[i];
    }
    return result;
}

double Vector::norm() const {
    return std::sqrt(dot(*this));
}

void Vector::print(std::string name) const {
    if (name != "") {
        std::cout << name << std::endl;
    }
    std::cout << "[ ";
    for (const auto& val : data) {
        std::cout << val << " ";
    }
    std::cout << "]" << std::endl;
}

} // namespace linalg
} // namespace FastFem