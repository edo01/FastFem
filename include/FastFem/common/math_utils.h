#pragma once
/******************************************************************************
 * Some good old friends
 *****************************************************************************/

#define POW2(x) ((x) * (x))
#define POW3(x) ((x) * (x) * (x))

#define MAX(A, B) ((B) > (A) ? (B) : (A))
#define MIN(A, B) ((B) < (A) ? (B) : (A))

#define CLAMP(x, low, high) \
	(((x) > (high)) ? (high) : (((x) < (low)) ? (low) : (x)))

#define PI 3.14159265358979323846
#define SQRT2 1.41421356237309504880
#define ONETHIRD 0.33333333333333333333
#define DEG2RAD 0.01745329251994329425

/* T must be floating point type */
template <typename T> inline bool approx_equal(T A, T B)
{
	/* All but two significand digits should match */
	T scaled_diff = (A - B) * 0.01f;
	return (A + scaled_diff == A || B - scaled_diff == B);
}

