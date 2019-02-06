#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ai.hh"
 
typedef struct { double x, y ,z; int group; } point_t, *point;
 
double randf(double m)
{
	return m * rand() / (RAND_MAX - 1.);
}
 
point gen_xy(int& count, double radius , std::vector<std::vector<double > >& data)
{
	double ang, r;
	double M_PI = 3.1415926535897932384626433;
	count = data.size();
	point p = static_cast<point_t*>(malloc(sizeof(point_t) * count) );
	point pt= static_cast<point_t*>(malloc(sizeof(point_t) * count) );
 
	std::size_t it;
	for (p = pt + count; p-- > pt;) {
		// ang = randf(2 * M_PI); // Генерация диска 
		// r = randf(radius);
		// p->x = r * cos(ang);
		// p->y = r * sin(ang);
		// p->z = r * sin(ang);
		p->x = data[it][0];
		p->y = data[it][1];
		p->z = data[it][2];
		
		++it;
	}
 
	return pt;
}
 
inline double dist2(point a, point b)
{
	double x = a->x - b->x, y = a->y - b->y , z = a->z - b->z;
	return x*x + y*y + z*z;
}
 //Возвращает кол-во ближайших
inline int nearest(point pt, point cent, int n_cluster, double *d2)
{
	int i, min_i;
	point c;
	double d, min_d;
	//for_n по кластерам
#	define for_n for(c = cent, i = 0; i < n_cluster; i++, c++)
	for_n {
		min_d = HUGE_VAL; // 
		min_i = pt->group;
		for_n {
			if (min_d > (d = dist2(c, pt))) {
				min_d = d; 
				min_i = i;
			}
		}
	}
	if (d2) *d2 = min_d;
	return min_i;
}
 // Нахождение центров 
void kpp(point pts, int len, point cent, int n_cent)
{
#	define for_len for (j = 0, p = pts; j < len; j++, p++)
	int i, j;
	int n_cluster;
	double sum, *d = static_cast<double*>(malloc(sizeof(double) * len) );
 
	point p, c;
	cent[0] = pts[ rand() % len ]; // первый выбирается случайным образом
	for (n_cluster = 1; n_cluster < n_cent; n_cluster++) {
		sum = 0;
		for_len {
			nearest(p, cent, n_cluster, d + j);
			sum += d[j];
		}
		sum = randf(sum);
		for_len {
			if ((sum -= d[j]) > 0) continue;
			cent[n_cluster] = pts[j];
			break;
		}
	}
	for_len p->group = nearest(p, cent, n_cluster, 0);
	free(d);
}
 
point lloyd(point pts, int len, int n_cluster)
{
	int i, j, min_i;
	int changed;
 
	point cent = static_cast<point_t*>(malloc(sizeof(point_t) * n_cluster) ),   p, c;
 
	/* assign init grouping randomly */
	//for_len p->group = j % n_cluster;
 
	/* or call k++ init */
	kpp(pts, len, cent, n_cluster);
	int v = 100;
	do {
		/* group element for centroids are used as counters */
		--v;
		
		for_n { 
		c->group = 0; 
		c->x = c->y =c->z= 0; 
		}
		for_len {
			c = cent + p->group;
			c->group++;
			c->x += p->x; 
			c->y += p->y;
			c->z += p->z;
		}
		for_n { 
		c->x /= c->group; 
		c->y /= c->group; 
		c->z /= c->group; 
		}
 
		changed = 0;
		/* find closest centroid of each point */
		for_len {
			min_i = nearest(p, cent, n_cluster, 0);
			if (min_i != p->group) {
				changed++;
				p->group = min_i;
				//std::cout<<"p->group = "<<p->group<<std::endl;
			}
		}
	//std::cout<< "	in 	"<<std::endl;
	} while(v > 0);
	//while (changed > (len >> 10)); /* stop when 99.9% of points are good */
	//while (changed > (len >> 5)); /* stop when 99.9% of points are good */
	//double a = len>>10;
	//std::cout<<"len >>10 = "<< a<< std::endl;
	for_n { c->group = i; }
 
	return cent;
}
 
void print_eps(point pts, int len, point cent, int n_cluster)
{
#	define W 400
#	define H 400
#	define L 400
	int i, j;
	point p, c;
	double min_x, max_x, min_y, max_y, min_z, max_z,scale, cx, cy, cz;
	double *colors = static_cast<double*>(malloc(sizeof(double) * n_cluster * 3) );
 
	for_n {
		colors[3*i + 0] = (3 * (i + 1) % 11)/11.;
		colors[3*i + 1] = (7 * i % 11)/11.;
		colors[3*i + 2] = (9 * i % 11)/11.;
	}
 
	max_x = max_y =max_z= -(min_x = min_y =min_z= HUGE_VAL);
	for_len {
		if (max_x < p->x) max_x = p->x;
		if (min_x > p->x) min_x = p->x;
		if (max_y < p->y) max_y = p->y;
		if (min_y > p->y) min_y = p->y;
		if (max_z < p->z) max_z = p->z;
		if (min_z > p->z) min_z = p->z;
	}
	scale = W / (max_x - min_x);
	if (scale > H / (max_y - min_y)) scale = H / (max_y - min_y);
	if (scale > L / (max_z - min_z)) scale = L / (max_z - min_z);
	cx = (max_x + min_x) / 2;
	cy = (max_y + min_y) / 2;
	cz = (max_z + min_z) / 2;
 
	printf("%%!PS-Adobe-3.0\n%%%%BoundingBox: -5 -5 %d %d\n", W + 10, H + 10);
	printf( "/l {rlineto} def /m {rmoveto} def\n"
		"/c { .25 sub exch .25 sub exch .5 0 360 arc fill } def\n"
		"/s { moveto -2 0 m 2 2 l 2 -2 l -2 -2 l closepath "
		"	gsave 1 setgray fill grestore gsave 3 setlinewidth"
		" 1 setgray stroke grestore 0 setgray stroke }def\n"
	);
	for_n {
		printf("%g %g %g setrgbcolor\n",
			colors[3*i], colors[3*i + 1], colors[3*i + 2]);
		for_len {
			if (p->group != i) continue;
			printf("%.3f %.3f c\n",
				(p->x - cx) * scale + W / 2,
				(p->y - cy) * scale + H / 2,
				(p->z - cz) * scale + L /2
				);
		}
		printf("\n0 setgray %g %g s\n",
			(c->x - cx) * scale + W / 2,
			(c->y - cy) * scale + H / 2,
			(c->z - cz) * scale + L /2
			);
	}
	printf("\n%%%%EOF");
	free(colors);

}
 
 
 void print_xyz(point pts, int len, point cent, int n_cluster , std::vector<std::vector<double > >& part , std::vector<std::vector<double > >& center, 
					std::vector<std::vector<std::vector<double> > >& cluster)
	{

		int i, j;
		point p, c;
		cluster.resize(n_cluster);
		for_n {
			for_len {
				//if (p->group != i) continue;
				// printf("%.3f %.3f c\n",
					// (p->x - cx) * scale + W / 2,
					// (p->y - cy) * scale + H / 2,
					// (p->z - cz) * scale + L /2
					// );
					if (p->group == i){
							cluster[i].push_back(std::vector<double> {p->x , p->y, p->z} );
					}
					part.push_back(std::vector<double> {p->x , p->y, p->z} );
			}
			// printf("\n0 setgray %g %g s\n",
				// (c->x - cx) * scale + W / 2,
				// (c->y - cy) * scale + H / 2,
				// (c->z - cz) * scale + L /2
				// );
				
				center.push_back(std::vector<double> {c->x , c->y, c->z} );
				
		}
	#undef for_n
	#undef for_len
	
	
	}

	
	double Average_Vector(std::vector<double>& vec){
		double sum = 0.;
		for(std::size_t i =0; i < vec.size() ; ++i){
			sum +=vec[i];
		}
		return sum / vec.size();
	}

 
 void RemoveRow(std::vector<std::vector<double> >& Matrix, std::size_t row){
	 
	 std::size_t size = Matrix.size();
	 //std::cout<<size<<std::endl;
	 for(std::size_t i = 0; i < size ; ++i){
		
			 if(i == row){
				 Matrix[i][0] = 0.;//Matrix[i+1][0];
				 Matrix[i][1] = 0.;//Matrix[i+1][1];
				 Matrix[i][2] = 0.;// Matrix[i+1][2];
		
			  }
	 }
	 //Matrix.resize(size - 1);
 }
 
 void Resizemartix(std::vector<std::vector<double> >& Matrix){
	 for(std::size_t i = 0 ; i < Matrix.size(); ++i){
		 if(Matrix[i][0] == 0. && Matrix[i][1] == 0. && Matrix[i][0] == 0.){
			 for(std::size_t j = i; j < Matrix.size()-1 ; ++j){
				 Matrix[j][0] = Matrix[j+1][0];
				 Matrix[j][1] = Matrix[j+1][1];
				 Matrix[j][2] = Matrix[j+1][2];
			 }
			 Matrix.resize(Matrix.size() - 1);
			 --i;
		 }
	 }
 }
 
 void cluster_MinDist(std::vector<std::vector<std::vector<double> > >& cluster, std::vector<std::vector<double > >& center){
	
	double dx , dy , dz;
	// std::cout<<"prev center"<<std::endl;
	// ai::printMatrix(center);
	for (std::size_t i = 0 ; i < cluster.size(); ++i){
		std::vector<double> dist;
		for(std::size_t j = 0 ; j < cluster[i].size() ; ++j){
			for(std::size_t k = j+1; k < cluster[i].size(); ++k ){
				dx = cluster[i][j][0] - cluster[i][k][0];
				dy = cluster[i][j][1] - cluster[i][k][1];
				dz = cluster[i][j][2] - cluster[i][k][2];
				
				dist.push_back( std::sqrt(dx*dx+dy*dy+dz*dz) );
				//std::cout<<"."<<std::endl;
			}
		}
		//ai::printVector(dist);
		//break;
		//std::cout<<Average_Vector(dist)<<std::endl;
		std::cout<<ai::min(dist)<<std::endl;
		 //if(Average_Vector(dist) > 45.){
		if(ai::min(dist) > 1.5){
			
			RemoveRow( center, i );
			//--i;
		}
			
	}
		 Resizemartix(center);
		// std::cout<<"resize center"<<std::endl;
		// ai::printMatrix(center);
 
 }
 
 