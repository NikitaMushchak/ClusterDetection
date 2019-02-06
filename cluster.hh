#pragma once
#include <iostream>
#include <vector>
#include <sstream>
#include <map>
#include <cmath>
#include <cfloat>
#include <string>
#include <string.h>
#include <math.h>
#include <float.h>
#include <ctype.h>
#include <utility>
#include "vect3d.h"






inline bool ut_isnan(double x)
{
  return isnan(x);
}

inline int ut_isinf(double x)
{
  return !isfinite(x);
}

//-----------------------------------------------------------------------------------------//
// Структура, описывающая частицу: координата, скорость, сила, масса частицы.
struct TParticle 
{
  Vect3D r, v, f;
  double m;
  int nint;
  double Q, T;
  //TCell_info base_cell_info;
  TParticle() : r(VECT3D_0), v(VECT3D_0), f(VECT3D_0), m(1.0), Q(0.0), T(0.0) {}
  TParticle(const Vect3D& zr) : r(zr), v(VECT3D_0), f(VECT3D_0), m(1.0), Q(0.0), T(0.0) {}
  //TParticle(const TShort_particle &p);
  //TParticle(const TShort_quasi &p);
};		 

struct TDump_config_info
{
  TDump_config_info() :
    iter(0),
    start_time(0),
    save_time(0),
    clusters_save_time(0),
    start_span_radius(0)
    {}

  int     iter;
  double  start_time; 
  double  save_time;
  double  clusters_save_time; 
  double  start_span_radius;
};

inline void swap_double(double &x)
{
  char *p = reinterpret_cast<char*>(&x);
  char q[8] = {p[7], p[6], p[5], p[4], p[3], p[2], p[1], p[0]};
  memcpy((void*)(&x), (void*)q, 8); 
}





int load_part_dump(std::string fname, 
                   std::vector<TParticle> &particles, 
                   TDump_config_info & info,
                   const int swap_bytes
                 );
				 
				 
				 

//#endif		





