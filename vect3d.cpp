//-----------------------------------------------------------------------------------
//
//	3D Vector class
//
//	Anton M. Krivtsov
//
//	11.04.2001
//
//-----------------------------------------------------------------------------------

#include "vect3d.h"
//#include "util.h"


inline double rand(const double min, const double max)	//	Random double value from min to max
{ 
	return min + (max - min) * (double)rand() / RAND_MAX; 
}



//-----------------------------------------------------------------------------------

double Vect3D::Abs() const	{ return sqrt(Sqr()); }

//-----------------------------------------------------------------------------------

void Vect3D::Swap(int mode)				//	 mode = 0,1,2;  default mode = 1
{
	if (!mode) return;
	double zz = z;
	if (mode == 1)	{ z = x; x = y; y = zz; return; } 
	if (mode == 2)	{ z = y; y = x; x = zz;	return; }
}

Vect3D Vect3D::GetSwapped(int mode) const	//	 mode = 0,1,2;  default mode = 1
{
	Vect3D r = *this; r.Swap(mode); return r;
}

//-----------------------------------------------------------------------------------

void Vect3D::SetRand(const double x_min, const double y_min, const double z_min, const double x_max, const double y_max, const double z_max) 
{
	x = rand(x_min, x_max);
	y = rand(y_min, y_max);
	z = rand(z_min, z_max);
}

//-----------------------------------------------------------------------------------

void Vect3D::SetRand(const double r_max) 
{
	double rr_max = r_max * r_max;
	do SetRand(-r_max, -r_max,0 , r_max,r_max , 0);
	while (Sqr() > rr_max);
}

//-----------------------------------------------------------------------------------
//	private
//-----------------------------------------------------------------------------------

void Vect3D::Rotate(double &x, double &y, const double sin_f, const double cos_f)
{ 
	double x1;
	x1 =  x * cos_f - y * sin_f; 
	y  =  x * sin_f + y * cos_f; 
	x  =  x1;
}

//-----------------------------------------------------------------------------------

Vect3D ConvertXY(double ro,double fi,double zz)
  {
    long double x,y,z;
    x=ro*cos(fi);
    y=ro*sin(fi);
    z=zz;
    Vect3D A(x,y,z);
    return A;
  }
  
//-----------------------------------------------------------------------------------
  
