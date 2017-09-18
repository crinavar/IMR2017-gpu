#include <iostream>
#include <sstream>
#include "c_mesh.h"
using namespace std;

struct timeval t_ini, t_fin;

inline void startTimer(){
    gettimeofday(&t_ini, NULL); //Tiempo de Inicio
}
inline double stopTimer(){
    gettimeofday(&t_fin, NULL); //Tiempo de Termino
    return (double)(t_fin.tv_sec + (double)t_fin.tv_usec/1000000) - (double)(t_ini.tv_sec + (double)t_ini.tv_usec/1000000);
}

int main(int argc, char **argv){
	double time=0;
	int ef=0;
	stringstream ss;
	if(argc < 2){
		printf("must run as ./app <offmesh_filename>\n");
		exit(1);
	}
	printf("loading mesh.....................");fflush(stdout);
	c_malla *m = new c_malla(argv[1]);
	printf("done\n"); fflush(stdout);
	//printf("min angle......."); fflush(stdout);
	double mangle = m->min_angle();
	//printf("%e rad (%f degrees)\n", mangle, 180.0*mangle/LAWSON_PI);
	printf("transforming into Delaunay.......");fflush(stdout);

	startTimer();
		ef = m->cpu_improve();
	time = stopTimer();
	printf("done: %.5g secs\n", time);
	printf("%f\n", time);

	ss << argv[1] << ".law.off";
	//printf("saving as [%s].......", ss.str().c_str()); fflush(stdout);	
	//m->export_off(ss.str().c_str());
	//printf("ok\n");fflush(stdout);
	//printf("appending result to [lawson_random.dat].......");fflush(stdout);
	//m->append_result(m->getNumVertices(), time, ef, "lawson_random.dat");
	//printf("ok\n");fflush(stdout);
	return 0;
	
}
