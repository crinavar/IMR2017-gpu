#include <stdio.h>
#include <cleap.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>
 
double delta_time(struct timeval t0, struct timeval t1);
 
int main(int argc, char* argv[]){
        cleap_mesh* m;
        struct timeval t0, t1;
        char buf[512];
        strcpy(buf, argv[1]);
        strcat(buf, ".mdt.off");
        cleap_init_no_render();
        m = cleap_load_mesh(argv[1]); // pass the mesh file, returns the mesh
        gettimeofday(&t0,NULL);
        cleap_delaunay_transformation(m, CLEAP_MODE_2D);
        gettimeofday(&t1,NULL);
        double dt = delta_time(t0, t1);
        printf("%f\n", dt);
    
	    cleap_save_mesh(m, buf);
	    cleap_end();
        return;
}

double delta_time(struct timeval t0, struct timeval t1){
	return (t1.tv_sec-t0.tv_sec) + ((double)(t1.tv_usec-t0.tv_usec)/1000000.0);
}

