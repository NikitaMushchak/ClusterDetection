#include <vector>
#include <cmath>
#include <cassert>
#include "cluster.hh"
#include "datareader.hh"
//#include "params.h"

const int header_length = 1024;


int is_inverse_byte_order()
{
  short x = 1 ;
  char *p = reinterpret_cast<char*> (&x);
return 0 != (int)(*p);}



int load_part_dump(std::string fname, 
                   std::vector<TParticle> &particles, 
                   TDump_config_info & info,
                   const int swap_bytes
                 )
{

  data_reader header;
  if(!header.open_file(fname))  return 0;
  int size = 0;
  header.get_int("size", size);
  header.get_int("dump_iter", info.iter);
  header.get_double("t_calc", info.start_time);
  header.get_double("t_save", info.save_time);
  header.get_double("t_clusters_save", info.clusters_save_time);
  header.get_double("start_span_radius", info.start_span_radius);

  particles.resize((int)size);
  FILE* f = fopen(fname.c_str(), "rb");
  
    char buf[header_length];
    fread(buf, header_length, 1, f);
  

  for (int i = 0; i < size; i++)
  {
    fread(&(particles[i].r), sizeof(double), 3, f);
    fread(&(particles[i].v), sizeof(double), 3, f);
    fread(&(particles[i].T), sizeof(double), 1, f);
    if (swap_bytes ^ !(is_inverse_byte_order())) 
    {
      swap_double(particles[i].r.x);
      swap_double(particles[i].r.y);
      swap_double(particles[i].r.z);
      swap_double(particles[i].v.x);
      swap_double(particles[i].v.y);
      swap_double(particles[i].v.z);
      swap_double(particles[i].T);
    }
    const Vect3D *r = &(particles[i].r);
    const Vect3D *v = &(particles[i].v);
    const double *t = &(particles[i].T);

      if(ut_isnan(r->x) || ut_isinf(r->x) ||
         ut_isnan(r->y) || ut_isinf(r->y) ||
         ut_isnan(r->z) || ut_isinf(r->z) ||
         ut_isnan(v->x) || ut_isinf(v->x) ||
         ut_isnan(v->y) || ut_isinf(v->y) ||
         ut_isnan(v->z) || ut_isinf(v->z) ||
          ut_isnan(*t) || ut_isinf(*t) ||
       r->x > 1e+10 || r->x < -1e+10 ||
       r->y > 1e+10 || r->y < -1e+10 ||
       r->z > 1e+10 || r->z < -1e+10 ||
       v->x > 1e+10 || v->x < -1e+10 ||
       v->y > 1e+10 || v->y < -1e+10 ||
       v->z > 1e+10 || v->z < -1e+10 ||
       (*t) > 1e+10 || (*t) < -1e+10 )

	   
      i--;
      size--;
      particles.pop_back();

      std::stringstream i_sstr;
      i_sstr << "load_part_dump :: wrong_particles_coords " << " p.x = " << r->x << " p.y = " << r->y << " p.z = " << r->z << " particles.size() = " << (int)particles.size();
	   
      //assert(0==1, i_sstr.str().c_str())
    

  }
  // char buff[512];

  // int nbufread = (int)fread(buff, sizeof(char), params.auxilinnary.end_part_dump_signature.length(), f);
  // buff[nbufread] = 0;

   fclose(f);

  // int cmp_res = strcmp(buff, params.auxilinnary.end_part_dump_signature.c_str());

  // if(cmp_res != 0)
  // {
    // string message = ((string)"Corrupted part dump file : " + fname);
    // //Log.PushMessage(message);
	// std::cout<<message<<std::endl;
    // //assert(cmp_res == 0, message.c_str());
  // }
	return 1;
}



 


  
  
  
  
  
