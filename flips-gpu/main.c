#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>
#include <cleap.h>
 
double delta_time(struct timeval t0, struct timeval t1);
 
int main(int argc, char* argv[]){
        if(argc != 2){
            fprintf(stderr, "run as ./mdt meshfile\n");
            exit(EXIT_FAILURE);
        }

        cleap_mesh* m;
        struct timeval t0, t1;
        char buf[512];
        strcpy(buf, argv[1]);
        strcat(buf, ".mdt.off");
        cleap_init_no_render();
        printf("loading mesh '%s'......", argv[1]); fflush(stdout);
        m = cleap_load_mesh(argv[1]); // pass the mesh file, returns the mesh
        printf("done\n"); fflush(stdout);
        gettimeofday(&t0,NULL);

        printf("delaunay transformation......", argv[1]); fflush(stdout);
        cleap_delaunay_transformation(m, CLEAP_MODE_2D);
        printf("done\n"); fflush(stdout);
        gettimeofday(&t1,NULL);
        double dt = delta_time(t0, t1);
        printf("%f\n", dt);
    
	    //cleap_save_mesh(m, buf);
	    cleap_end();
        return 0;
}

double delta_time(struct timeval t0, struct timeval t1){
	return (t1.tv_sec-t0.tv_sec) + ((double)(t1.tv_usec-t0.tv_usec)/1000000.0);
}

