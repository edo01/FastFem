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

Vector Vector::axpby(double a, const Vector& x, double b, const Vector& y) const {
    Vector result(data.size());
    #pragma omp parallel for
    for (std::size_t i = 0; i < data.size(); ++i) {
        result.data[i] = a * x.data[i] + b * y.data[i];
    }
    return result;
}

void Vector::axpby(double a, const Vector& x, double b, Vector& y) const {
    #pragma omp parallel for
    for (std::size_t i = 0; i < data.size(); ++i) {
        y.data[i] = a * x.data[i] + b * y.data[i];
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

void Vector::print() const {
    std::cout << "[ ";
    for (const auto& val : data) {
        std::cout << val << " ";
    }
    std::cout << "]" << std::endl;
}

} // namespace linalg
} // namespace FastFem