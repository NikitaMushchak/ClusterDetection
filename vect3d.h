#pragma once

#ifndef ___Vect3D_h___
#define ___Vect3D_h___

#include <math.h>
//#include "util.h"


#define Vect3D_ const Vect3D&

class Box3D;

//-----------------------------------------------------------------------------------

class Vect3D
{
public:
	struct	 SFLOAT								{ float x;	float y;  float z; };
	struct	 SDOUBLE							{ double x;	double y; double z; };

	Vect3D()									{x = 0; y = 0; z = 0;}				
	Vect3D( const double  x,  const double  y,  const double  z)					{ Set(x, y, z); } 
	Vect3D(const SFLOAT&  r)					{ x = r.x; y = r.y; z = r.z; }
	Vect3D(const SDOUBLE& r)					{ x = r.x; y = r.y; z = r.z; }

	void Set( const double  x,  const double  y,  const double  z)				{ this->x = x; this->y = y; this->z = z; }

	void SetRand( const double  x_min,  const double  y_min,  const double  z_min,  const double  x_max,  const double  y_max,  const double  z_max);
	void SetRand(const Box3D&);					//	defined in Box3D.cpp
	void SetRand( const double  r_max);					

	operator SFLOAT()					const	{ SFLOAT r = {float(x), float(y), float(z)}; return r; }
	operator SDOUBLE()					const	{ SDOUBLE r = {x, y, z}; return r; }

	Vect3D& operator +=(Vect3D_ r)				{ x += r.x;	y += r.y; z += r.z; return *this; }
	Vect3D& operator -=(Vect3D_ r)				{ x -= r.x;	y -= r.y; z -= r.z; return *this; }
	Vect3D& operator *=(const double k)					{ x *= k;	y *= k;   z *= k;	return *this; }
  Vect3D& operator /=(const double k)					{ x /= k;	y /= k;   z /= k;	return *this; }

	Vect3D operator -()					const	{ Vect3D r(-x, -y, -z); return r; }

	int operator ==(Vect3D_ r)			const	{ return x==r.x && y==r.y && z==r.z; }
	int operator !=(Vect3D_ r)			const	{ return !operator==(r); }
	int IsNull()						const	{ return !(x || y || z); }

	void FlipX()								{ x = -x; }
	void FlipY()								{ y = -y; }
	void FlipZ()								{ z = -z; }

	void RotateZ( const double  angle)						{ Rotate(x, y, sin(angle), cos(angle)); }
	void RotateX( const double  angle)						{ Rotate(y, z, sin(angle), cos(angle)); }
	void RotateY( const double  angle)						{ Rotate(z, x, sin(angle), cos(angle)); }

	void RotateX( const double  sin_f,  const double  cos_f)			{ Rotate(y, z, sin_f, cos_f); }
	void RotateY( const double  sin_f,  const double  cos_f)			{ Rotate(z, x, sin_f, cos_f); }
	void RotateZ( const double  sin_f,  const double  cos_f)			{ Rotate(x, y, sin_f, cos_f); }

	void RotateZXZ( const double  phi,  const double  theta,  const double  psi)	{ RotateZ(phi); RotateX(theta); RotateZ(psi); }

	void   Swap		  (int mode = 1);			// mode = 0,1,2
	Vect3D GetSwapped (int mode = 1)	const;	// mode = 0,1,2

	double Sqr()							const	{ return x * x + y * y + z * z; }
	double Abs()							const;

	double x, y, z;

private:
	static void Rotate(double  &x, double  &y,  const double  sin_f,  const double  cos_f);
};

//-----------------------------------------------------------------------------------

static const Vect3D VECT3D_0(0, 0, 0);		// zero vector

//-----------------------------------------------------------------------------------

static inline Vect3D operator +(Vect3D_ r1, Vect3D_ r2)	{ Vect3D rr(r1);	return rr += r2;	}
static inline Vect3D operator -(Vect3D_ r1, Vect3D_ r2)	{ Vect3D rr(r1);	return rr -= r2;	}
static inline Vect3D operator *(Vect3D_ r,   const double  k)		{ Vect3D rr(r);		return rr *= k;		}
static inline Vect3D operator *( const double  k,		Vect3D_ r)	{					return r * k;		}
static inline Vect3D operator /(Vect3D_ r,   const double  k)		{ Vect3D rr(r);		return rr /= k;		}



static inline Vect3D operator %(Vect3D_ r1, Vect3D_ r2)	// vector product
{ return Vect3D(
	r1.y * r2.z - r1.z * r2.y,
	r1.z * r2.x - r1.x * r2.z,
	r1.x * r2.y - r1.y * r2.x);	
}

static inline double	 operator *(Vect3D_ r1, Vect3D_ r2)	{ return r1.x * r2.x + r1.y * r2.y + r1.z * r2.z; }

//-----------------------------------------------------------------------------------
Vect3D ConvertXY(double ro,double fi,double zz);
  
#endif //___Vect3D_h___

