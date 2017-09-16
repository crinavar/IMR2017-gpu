/**
 *@author Cristobal Navarro
 *@contact axischire@gmail.com
 *@detail Clase c_malla
 *@creation 1354 Jan 25 2010.
 */
#ifndef C_MESH_H
#define C_MESH_H

#define EPSILON_DELAUNAY 	0.000001f
#define LAWSON_PI		3.14159265f

#include <cstdio>
#include <string>
#include <set>
#include <map>
#include <vector>
#include <stack>
#include <limits>
#include <cstdlib>
#include <cmath>
#include <tr1/unordered_map>
#include "structs.h"

class c_malla{

    public:

        edges edge_data;
        mesh local_mesh;
        Vertex* vertexesMalla;
        unsigned int* triangulosMalla;
        unsigned int* triangle_edges;
        unsigned int* edge_marks;
        int* collisionsMap;
        int edgesProcesados;
        int paresSize;
	int numVertexes, numFaces, numEdges, numIndexes;

        std::stack<unsigned int> edge_stack;
        unsigned int* cpu_trirel;
        
        float minx, miny, minz, maxx, maxy, maxz;

//funciones
        c_malla();
        c_malla(char* filename);
        c_malla(int numVertexes, int numFaces, int numEdges);
        void setCuantities( int numVertexes, int numFaces, int numEdges );
        int getNumVertices();
        int getNumFaces();
        int getNumEdges();
        int getNumIndexes();
	double min_angle();
        void resetMinMax();
        void printTriangleArray();
        void export_off(const char *filename);
        void append_result(int n, double t, int ef, const char *filename);
        ~c_malla();

//algoritmos por version-secuencial
        int cpu_improve();
        void cpu_improve_3d();
        int cpu_fix_edge(int i, int rt);
        void cpu_fix_add_neighbors( int i );
        int cpu_2d_delaunay(int i);
        int cpu_3d_delaunay(int i, const float limit_angle );
        void cpu_edge_flip( int i);
        void make_stack();
        void print_edge(int i);

};

#endif
