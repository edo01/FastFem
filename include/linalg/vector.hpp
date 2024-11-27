#ifndef VECTOR_HPP
#define VECTOR_HPP

#include <cmath>

namespace linalg
{
struct Vector 
{
	double x;
	double y;
	double z;
};

double dot(struct Vector V, struct Vector W);

double norm(struct Vector V);

struct Vector cross(struct Vector V, struct Vector W);

}

#endif // VECTOR_HPP