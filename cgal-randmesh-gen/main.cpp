#include "saver.h"

int main (int argc, char **argv){   
	int N = 1000000; if( argc > 1 )N = atoi(argv[1]); // number of points
    if( argc != 3 ){
        fprintf(stderr, "error: bad arguments: try running this way --> ./rgen <number_of_points> <outfile>\n");
        exit(1);
    }
	Rnd_points_it   rand_it( 0.499 );
	Points_iterator begin(rand_it), end(rand_it,N);
	CT ct;
	CGAL::Timer cost; cost.reset();cost.start();
	std::cout << "constrained triangulation of " << N << " points in 2D unit square" << std::flush;
	ct.insert_constraint(Point(-0.5, -0.5), Point(0.5, -0.5));
	ct.insert_constraint(Point(0.5, -0.5), Point(0.5, 0.5));	
	ct.insert_constraint(Point(0.5, 0.5), Point(-0.5, 0.5));	
	ct.insert_constraint(Point(-0.5, 0.5), Point(-0.5, -0.5));	
	ct.insert(begin,end);cost.stop();
	std::cout << " done in "<<cost.time()/N*1000000<<" micro-seconds per point." << std::endl;
	cost.reset();cost.start();
	ct.is_valid();
	std::cout << " checked in "<<cost.time()/N*1000000<<" micro-seconds per point." << std::endl;
	std::map<CT::Vertex_handle, int> vmap;
	int index = 0;
	for (CT::Finite_vertices_iterator vit = ct.finite_vertices_begin(), end = ct.finite_vertices_end(); vit != end; ++vit)
		vmap[vit] = index++;
	save_mesh(ct, vmap, argv[2]);
	return 0;
}
