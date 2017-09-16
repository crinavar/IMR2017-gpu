#ifndef STRUCTS_H
#define STRUCTS_H
#include <sys/time.h>

typedef std::pair<int, int> pairArco;

typedef struct{
	float x, y, z, w;
} float4;

typedef struct{
	float x, y, z;
} float3;

typedef struct{
	int x, y;
} int2;

typedef struct{

    //cada float son 4 bytes
    float x, y, z;
    float nx, ny, nz;
    float r, g, b, a;
    //grow factor
    float gFactor;

} Vertex;

typedef struct{

    int id;
    int n1, n2;
    int a1, a2;
    int b1, b2;
    int op1, op2;
} Edge;


typedef struct{

	float4 *v;
	float4 *n;
	float4 *c;

} mesh;

typedef struct{
	int2 *n;
	int2 *a;
	int2 *b;
	int2 *op;
} edges;

#endif
