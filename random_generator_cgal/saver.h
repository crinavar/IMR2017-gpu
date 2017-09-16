#ifndef _SAVER_H_
#define _SAVER_H_

#include <iostream>
#include <fstream>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/iterator.h>
#include <CGAL/Timer.h>
#include <vector>
#include <map>
#include <sys/time.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel		K;
typedef CGAL::Triangulation_vertex_base_2<K>				Vb;
typedef CGAL::Constrained_triangulation_face_base_2<K>     		Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb>			TDS;
typedef CGAL::Exact_predicates_tag					Itag;
typedef CGAL::Constrained_triangulation_2<K, TDS, Itag>			CT;
typedef CT::Point	    Point;
typedef CGAL::Random_points_in_square_2<Point>          	Rnd_points_it;
typedef CGAL::Counting_iterator<Rnd_points_it, Point>   	Points_iterator;

double diff_time(timeval t0, timeval t1){
	return (t1.tv_sec-t0.tv_sec) + (double(t1.tv_usec-t0.tv_usec)/1000000.0);
}
void save_mesh(	CT &t,	std::map<CT::Vertex_handle, int> vmap, const char *filename){
	printf("saving constrained triangulation.......");fflush(stdout);
	std::ofstream fw(filename);
	fw << "OFF" << std::endl;
	fw << t.number_of_vertices() << " " << t.number_of_faces() << " " << 0 << std::endl;
	//for(CDT::Finite_vertices_iterator it = t.finite_vertices_begin(); it != t.finite_vertices_end(); it++)
	//	fw << *it << " " << 0.0 << std::endl;
	for (CT::Finite_vertices_iterator vit = t.finite_vertices_begin(), end = t.finite_vertices_end(); vit != end; ++vit)
		fw << *vit << " " << 0.0 << std::endl;
	for(CT::Finite_faces_iterator it = t.finite_faces_begin(); it != t.finite_faces_end(); it++)
		fw << 3 << " " << vmap.find(it->vertex(0))->second << " " << 
				  vmap.find(it->vertex(1))->second << " " << 
				  vmap.find(it->vertex(2))->second << std::endl;
	printf("ok\n");fflush(stdout);
	fw.close();
}
#endif
