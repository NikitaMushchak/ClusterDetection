#include <iostream>
#include <vector>
#include "ai.hh"
#include "cluster.hh"
#include "kalgorithm.h"



//#define PTS 100000
#define K 30
 //Число кластеров

/// PTS - кол-во частиц

int main(){
	//TParticle Tparticles;
	
	std::vector<TParticle> particles;
	
	std::string fname("./wr000_r000__ws008_s004_part_dump.dat");
	int swap_bytes = 0;
	TDump_config_info info;
	
	 load_part_dump(fname, 
                   particles, 
                    info,
                    swap_bytes);
	
	std::cout<< particles.size()<<std::endl;
	std::vector<std::vector<double> > data;
	
	std::vector <TParticle> ::iterator p;
	
			for (p = particles.begin(); p < particles.end(); p++) {
					
				data.push_back(std::vector<double> {p->r.x, p->r.y, p->r.z});
			}
			
			//ai::printMatrix(data);
			const double radius = 1;
			int PTS;
			//ai::saveA3R("./save",data,radius);
			std::vector<std::vector<double > > part;
			
			std::vector<std::vector<std::vector<double > > >cluster;
			
			int i;
			
			std::string num ;
			
			point v = gen_xy(PTS, 10, data);
			data.clear();
			
				
				std::vector<std::vector<double > > center;	
				point c = lloyd(v, PTS, K);
			
				print_xyz(v, PTS, c, K ,  part ,  center, cluster);
				
				num = std::to_string(i);
				
				ai::saveXYZ(num + "cent", center);
				
				cluster_MinDist( cluster, center);
				
				std::vector<std::vector<double> > cl;
				for(std::size_t i =0 ; i < cluster[1].size(); ++i)
					cl.push_back(std::vector<double>{cluster[1][i][0], cluster[1][i][1], cluster[1][i][2] } 	); 
					
				ai::saveXYZ( "clasterrrrr", cl);
				
				std::cout<<"cluster1 size = "<<cluster[1].size()<<std::endl;
					
				std::vector<std::vector<double> > ce;
				
				ce.push_back(std::vector<double>{center[1][0] , center[1][1]	,	center[1][2]	}	);
				
				ai::saveXYZ( "middle", ce);
				
				
				
				if(i==0)
				ai::saveXYZ( "part", part);
							
				ai::saveXYZ(num + "clean_cent", center);
				delete(c);
			
			
			std::cout<<"cluster size :"<<cluster.size();
		//	ai::saveXYZ("cluster", cluster);
			
	// free(v); free(c);
			
			
        
	return 1;
}