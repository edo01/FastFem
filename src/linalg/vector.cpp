# include "linalg/vector.hpp"

namespace linalg
{

double dot(struct Vector V, struct Vector W) 
{
	return V.x*W.x + V.y*W.y + V.z*W.z;
}

double norm(struct Vector V) 
{ 
	return sqrt(dot(V, V)); 
};

struct Vector cross(struct Vector V, struct Vector W)
{
	Vector res;
	// cross product between V and W
	res.x = V.y*W.z - V.z*W.y; 
	res.y = V.z*W.x - V.x*W.z;
	res.z = V.x*W.y - V.y*W.x;

	return res;
}

} // linalg