#include "c_mesh.h"

int2 make_int2(int a, int b){
	int2 v = {a, b};
	return v;
}

float3 make_float3(float a, float b, float c){
	float3 v = {a, b, c};
	return v;
}

float4 make_float4(float a, float b, float c, float d){
	float4 v = {a, b, c, d};
	return v;
}


c_malla::c_malla(){
}

c_malla::c_malla(char *filename){

    //Lineas necesarias para que scanf lea archivo en computadores seteados en otro lenguaje.
    setlocale(LC_NUMERIC, "POSIX");
    FILE *off = fopen(filename,"r");
    char line[255];
    int num_puntos, num_caras, num_arcos;

    fscanf(off,"%s\n",line);
    while(true){
        fgets(line,255,off);
        if (line[0]!='#')break;
    }
    sscanf(line,"%d %d %d\n",&num_puntos,&num_caras,&num_arcos);

    //c_malla *c_malla = 0;
    //c_malla = new c_malla( num_puntos, num_caras, num_arcos);
    this->setCuantities( num_puntos, num_caras, this->paresSize );
    this->resetMinMax();

    //metodo array of structs
    this->vertexesMalla = (Vertex*)malloc(sizeof(Vertex)*num_puntos);
    this->triangulosMalla = (unsigned int*)malloc(sizeof(unsigned int)*num_caras*3);

    //metodo struct of arrays
    this->local_mesh.v = (float4*)malloc(sizeof(float4)*num_puntos);
    this->local_mesh.n = (float4*)malloc(sizeof(float4)*num_puntos);
    this->local_mesh.c = (float4*)malloc(sizeof(float4)*num_puntos);

    //--PARSING VERTEXES
    float cont=1.0f;
    for(int i=0; i<num_puntos; i++) {

        fscanf(off,"%f %f %f\n",&this->vertexesMalla[i].x,&this->vertexesMalla[i].y,&this->vertexesMalla[i].z);
        this->vertexesMalla[i].nx = 0.0f;
        this->vertexesMalla[i].ny = 0.0f;
        this->vertexesMalla[i].nz = 0.0f;
        this->vertexesMalla[i].gFactor = 1.0f;

        //struct of arrays method
        //fscanf(off,"%f %f %f\n",&this->local_mesh.x[i],&this->local_mesh.y[i],&this->local_mesh.z[i]);

        this->local_mesh.v[i].x = this->vertexesMalla[i].x;
        this->local_mesh.v[i].y = this->vertexesMalla[i].y;
        this->local_mesh.v[i].z = this->vertexesMalla[i].z;
        this->local_mesh.v[i].w = 1.0f;

        this->local_mesh.n[i].x = 0.0f;
        this->local_mesh.n[i].y = 0.0f;
        this->local_mesh.n[i].z = 0.0f;
        this->local_mesh.n[i].w = 1.0f;

        if (this->vertexesMalla[i].x>this->maxx) this->maxx=this->vertexesMalla[i].x;
        if (this->vertexesMalla[i].y>this->maxy) this->maxy=this->vertexesMalla[i].y;
        if (this->vertexesMalla[i].z>this->maxz) this->maxz=this->vertexesMalla[i].z;
        if (this->vertexesMalla[i].x<this->minx) this->minx=this->vertexesMalla[i].x;
        if (this->vertexesMalla[i].y<this->miny) this->miny=this->vertexesMalla[i].y;
        if (this->vertexesMalla[i].z<this->minz) this->minz=this->vertexesMalla[i].z;

        if (this->vertexesMalla[i].x>this->maxx) this->maxx=this->local_mesh.v[i].x;
        if (this->vertexesMalla[i].y>this->maxy) this->maxy=this->local_mesh.v[i].y;
        if (this->vertexesMalla[i].z>this->maxz) this->maxz=this->local_mesh.v[i].z;
        if (this->vertexesMalla[i].x<this->minx) this->minx=this->local_mesh.v[i].x;
        if (this->vertexesMalla[i].y<this->miny) this->miny=this->local_mesh.v[i].y;
        if (this->vertexesMalla[i].z<this->minz) this->minz=this->local_mesh.v[i].z;

    }
    //-------------------
    int tipo_cara;
    float r,g,b;
    char resto[256];
    int ch;
    int face = 3;
    float3 normal;
    float3 v1,v2;
    cont=1.0f;
    //PARSING FACES & NORMALS
    std::map<pairArco, Edge> mapaEdges;
    std::map<pairArco, Edge>::iterator it;
    std::map<pairArco, Edge>::iterator it2;
    std::vector<Edge*> vectorIndexPair;
    unsigned int *triangle_edges = (unsigned int*)malloc(sizeof(unsigned int)*num_caras*3);
    int tri_edge_count = 0;

    Edge* auxEdgePointer;
    //progressDialog->set_title("Loading mesh...");
    for(int i=0; i<num_caras; i++) {
        fscanf(off,"%d",&tipo_cara);
        if( tipo_cara == 3 ){

            fscanf(off,"%d %d %d",&this->triangulosMalla[i*tipo_cara], &this->triangulosMalla[i*tipo_cara+1], &this->triangulosMalla[i*tipo_cara+2]);
            //Construyendo Aristas
            int j=0, k=1, op=2;
            it = mapaEdges.find(pairArco(this->triangulosMalla[i*tipo_cara+j], this->triangulosMalla[i*tipo_cara+k]));
            if( it != mapaEdges.end() ){
                auxEdgePointer = &it->second;
                auxEdgePointer->b1 = i*tipo_cara+j;
                auxEdgePointer->b2 = i*tipo_cara+k;
                auxEdgePointer->op2 = i*tipo_cara+op;
                triangle_edges[tri_edge_count] = auxEdgePointer->id;
            }
            else{
                it = mapaEdges.find(pairArco(this->triangulosMalla[i*tipo_cara+k], this->triangulosMalla[i*tipo_cara+j]));
                if( it != mapaEdges.end() ){
                    auxEdgePointer = &it->second;
                    auxEdgePointer->b1 = i*tipo_cara+k;
                    auxEdgePointer->b2 = i*tipo_cara+j;
                    auxEdgePointer->op2 = i*tipo_cara+op;
                    triangle_edges[tri_edge_count] = auxEdgePointer->id;
                }
                else{
                    auxEdgePointer = &mapaEdges[pairArco(this->triangulosMalla[i*tipo_cara+j], this->triangulosMalla[i*tipo_cara+k])];
                    auxEdgePointer->n1 = this->triangulosMalla[i*tipo_cara+j];
                    auxEdgePointer->n2 = this->triangulosMalla[i*tipo_cara+k];
                    auxEdgePointer->a1 = i*tipo_cara+j;
                    auxEdgePointer->a2 = i*tipo_cara+k;
                    auxEdgePointer->b1 = -1;
                    auxEdgePointer->b2 = -1;
                    auxEdgePointer->op1 = i*tipo_cara+op;
                    auxEdgePointer->op2 = -1;
                    auxEdgePointer->id = vectorIndexPair.size();
                    // agregar el id del edge a un arreglo de edges por triangulo.
                    triangle_edges[tri_edge_count] = auxEdgePointer->id;
                    vectorIndexPair.push_back( auxEdgePointer );
                }
            }
            j=0, k=2, op=1, tri_edge_count++;
            it = mapaEdges.find(pairArco(this->triangulosMalla[i*tipo_cara+j], this->triangulosMalla[i*tipo_cara+k]));
            if( it != mapaEdges.end() ){
                auxEdgePointer = &it->second;
                auxEdgePointer->b1 = i*tipo_cara+j;
                auxEdgePointer->b2 = i*tipo_cara+k;
                auxEdgePointer->op2 = i*tipo_cara+op;
                triangle_edges[tri_edge_count] = auxEdgePointer->id;
            }
            else{
                it = mapaEdges.find(pairArco(this->triangulosMalla[i*tipo_cara+k], this->triangulosMalla[i*tipo_cara+j]));
                if( it != mapaEdges.end() ){
                    auxEdgePointer = &it->second;
                    auxEdgePointer->b1 = i*tipo_cara+k;
                    auxEdgePointer->b2 = i*tipo_cara+j;
                    auxEdgePointer->op2 = i*tipo_cara+op;
                    triangle_edges[tri_edge_count] = auxEdgePointer->id;
                }
                else{
                    auxEdgePointer = &mapaEdges[pairArco(this->triangulosMalla[i*tipo_cara+j], this->triangulosMalla[i*tipo_cara+k])];
                    auxEdgePointer->n1 = this->triangulosMalla[i*tipo_cara+j];
                    auxEdgePointer->n2 = this->triangulosMalla[i*tipo_cara+k];
                    auxEdgePointer->a1 = i*tipo_cara+j;
                    auxEdgePointer->a2 = i*tipo_cara+k;
                    auxEdgePointer->b1 = -1;
                    auxEdgePointer->b2 = -1;
                    auxEdgePointer->op1 = i*tipo_cara+op;
                    auxEdgePointer->op2 = -1;
                    auxEdgePointer->id = vectorIndexPair.size();
                    triangle_edges[tri_edge_count] = auxEdgePointer->id;
                    vectorIndexPair.push_back( auxEdgePointer );
                }
            }
            j=1, k=2, op=0, tri_edge_count++;
            it = mapaEdges.find(pairArco(this->triangulosMalla[i*tipo_cara+j], this->triangulosMalla[i*tipo_cara+k]));
            if( it != mapaEdges.end() ){
                auxEdgePointer = &it->second;
                auxEdgePointer->b1 = i*tipo_cara+j;
                auxEdgePointer->b2 = i*tipo_cara+k;
                auxEdgePointer->op2 = i*tipo_cara+op;
                triangle_edges[tri_edge_count] = auxEdgePointer->id;
            }
            else{
                it = mapaEdges.find(pairArco(this->triangulosMalla[i*tipo_cara+k], this->triangulosMalla[i*tipo_cara+j]));
                if( it != mapaEdges.end() ){
                    auxEdgePointer = &it->second;
                    auxEdgePointer->b1 = i*tipo_cara+k;
                    auxEdgePointer->b2 = i*tipo_cara+j;
                    auxEdgePointer->op2 = i*tipo_cara+op;
                    triangle_edges[tri_edge_count] = auxEdgePointer->id;
                }
                else{
                    auxEdgePointer = &mapaEdges[pairArco(this->triangulosMalla[i*tipo_cara+j], this->triangulosMalla[i*tipo_cara+k])];
                    auxEdgePointer->n1 = this->triangulosMalla[i*tipo_cara+j];
                    auxEdgePointer->n2 = this->triangulosMalla[i*tipo_cara+k];
                    auxEdgePointer->a1 = i*tipo_cara+j;
                    auxEdgePointer->a2 = i*tipo_cara+k;
                    auxEdgePointer->b1 = -1;
                    auxEdgePointer->b2 = -1;
                    auxEdgePointer->op1 = i*tipo_cara+op;
                    auxEdgePointer->op2 = -1;
                    auxEdgePointer->id = vectorIndexPair.size();
                    vectorIndexPair.push_back( auxEdgePointer );
                    triangle_edges[tri_edge_count] = auxEdgePointer->id;
                }
            }
            tri_edge_count++;
        }
        else{
            //escribirEnLog("ABORTADO: la malla debe ser triangular. ");
            //escribirStatusLog(LOG_APP,  "[CUDA-MODE] ABORT: la malla debe ser triangular");
            free(this->local_mesh.v);
            free(this->local_mesh.n);
            free(this->local_mesh.c);
            free(this->triangulosMalla);
            fclose(off);
            printf("abort:: mesh must have only triangles");
            return;
        }
        Vertex p1 = this->vertexesMalla[this->triangulosMalla[i*face]];
        Vertex p2 = this->vertexesMalla[this->triangulosMalla[i*face+1]];
        Vertex p3 = this->vertexesMalla[this->triangulosMalla[i*face+2]];
        v1 = make_float3( p2.x - p1.x, p2.y - p1.y, p2.z - p1.z);
        v2 = make_float3( p3.x - p1.x, p3.y - p1.y, p3.z - p1.z);
        normal.x =   (v1.y * v2.z) - (v2.y * v1.z);
        normal.y = -((v1.x * v2.z) - (v2.x * v1.z));
        normal.z =   (v1.x * v2.y) - (v2.x * v1.y);
        //printf("   Normal= (%f, %f, %f)\n", normal.x, normal.y, normal.z);

        //printf("Seteando Normales en Nodo %i   (%f, %f, %f)", triangulosMalla[i*face], this->vertexesMalla[triangulosMalla[i*face]].x, this->vertexesMalla[triangulosMalla[i*face]].y, this->vertexesMalla[triangulosMalla[i*face]].z);
        this->vertexesMalla[this->triangulosMalla[i*face]].nx += normal.x;
        this->vertexesMalla[this->triangulosMalla[i*face]].ny += normal.y;
        this->vertexesMalla[this->triangulosMalla[i*face]].nz += normal.z;
        //printf("  N(%f, %f, %f)\n", this->vertexesMalla[triangulosMalla[i*face]].nx, this->vertexesMalla[triangulosMalla[i*face]].ny, this->vertexesMalla[triangulosMalla[i*face]].nz );

        //printf("Seteando Normales en Nodo %i   (%f, %f, %f)", triangulosMalla[i*face+1], this->vertexesMalla[triangulosMalla[i*face+1]].x, this->vertexesMalla[triangulosMalla[i*face+1]].y, this->vertexesMalla[triangulosMalla[i*face+1]].z);
        this->vertexesMalla[this->triangulosMalla[i*face+1]].nx += normal.x;
        this->vertexesMalla[this->triangulosMalla[i*face+1]].ny += normal.y;
        this->vertexesMalla[this->triangulosMalla[i*face+1]].nz += normal.z;
        //printf("  N(%f, %f, %f)\n", this->vertexesMalla[triangulosMalla[i*face+1]].nx, this->vertexesMalla[triangulosMalla[i*face+1]].ny, this->vertexesMalla[triangulosMalla[i*face+1]].nz );

        //printf("Seteando Normales en Nodo %i   (%f, %f, %f)", triangulosMalla[i*face+2], this->vertexesMalla[triangulosMalla[i*face+2]].x, this->vertexesMalla[triangulosMalla[i*face+2]].y, this->vertexesMalla[triangulosMalla[i*face+2]].z);
        this->vertexesMalla[this->triangulosMalla[i*face+2]].nx += normal.x;
        this->vertexesMalla[this->triangulosMalla[i*face+2]].ny += normal.y;
        this->vertexesMalla[this->triangulosMalla[i*face+2]].nz += normal.z;
        //printf("  N(%f, %f, %f)\n", this->vertexesMalla[triangulosMalla[i*face+2]].nx, this->vertexesMalla[triangulosMalla[i*face+2]].ny, this->vertexesMalla[triangulosMalla[i*face+2]].nz );

        //!METODO ARRAY OF STRUCTS
        this->local_mesh.n[this->triangulosMalla[i*face]].x += normal.x;
        this->local_mesh.n[this->triangulosMalla[i*face]].y += normal.y;
        this->local_mesh.n[this->triangulosMalla[i*face]].z += normal.z;

        this->local_mesh.n[this->triangulosMalla[i*face+1]].x += normal.x;
        this->local_mesh.n[this->triangulosMalla[i*face+1]].y += normal.y;
        this->local_mesh.n[this->triangulosMalla[i*face+1]].z += normal.z;

        this->local_mesh.n[this->triangulosMalla[i*face+2]].x += normal.x;
        this->local_mesh.n[this->triangulosMalla[i*face+2]].y += normal.y;
        this->local_mesh.n[this->triangulosMalla[i*face+2]].z += normal.z;

    }
    this->edgesProcesados = 0;
    this->paresSize = vectorIndexPair.size();
    this->edge_data.n = (int2*)malloc( sizeof(int2)*this->paresSize );
    this->edge_data.a = (int2*)malloc( sizeof(int2)*this->paresSize );
    this->edge_data.b = (int2*)malloc( sizeof(int2)*this->paresSize );
    this->edge_data.op = (int2*)malloc( sizeof(int2)*this->paresSize );
    //set array of edge_marks
    this->edge_marks = (unsigned int*)malloc( sizeof(unsigned int)*this->paresSize );
    //set triangle edges for secuential algorithm
    this->triangle_edges = triangle_edges;
    for( int i=0; i<this->paresSize; i++ ){
            this->edge_data.n[i] = make_int2(vectorIndexPair[i][0].n1, vectorIndexPair[i][0].n2);
            this->edge_data.a[i] = make_int2(vectorIndexPair[i][0].a1, vectorIndexPair[i][0].a2);
            this->edge_data.b[i] = make_int2(vectorIndexPair[i][0].b1, vectorIndexPair[i][0].b2);
            this->edge_data.op[i] = make_int2(vectorIndexPair[i][0].op1, vectorIndexPair[i][0].op2);
            this->edge_marks[i] = 0;
    }
    //printf("triangle edges success!!!!\n");
/*
    for( int i=0; i<this->paresSize; i++ ){
            printf("edge %i:\n", i);
            printf("n = (%i, %i)\t", this->edge_data.n[i].x, this->edge_data.n[i].y);
            printf("a = (%i, %i)\t", this->edge_data.a[i].x, this->edge_data.a[i].y);
            printf("b = (%i, %i)\n--------\n", this->edge_data.b[i].x, this->edge_data.b[i].y);
    }
*/
    //printf("TIEMPO caras: %.5g[ms]\n", stopTimer()*1000.0 );

    //----OPERAR CON LOS DATOS
    //----
	fclose(off);
	setlocale(LC_NUMERIC, "");
    //printf("[v=%i, e=%i, f=%i].......", num_puntos, this->paresSize, num_caras);
    this->setCuantities( num_puntos, num_caras, this->paresSize );
/*
    printf("Mesh= v=%i\n", num_puntos);
    for(int i=0; i<num_puntos; i++){
        printf("v(x,y,z,w) = (%f, %f, %f, %f)\n", this->local_mesh.v[i].x, this->local_mesh.v[i].y, this->local_mesh.v[i].z, this->local_mesh.v[i].w);
    }
*/
    mapaEdges.clear();
    vectorIndexPair.clear();
    this->make_stack(); //creando stack de edges a procesar
}


c_malla::c_malla(int numVertexes, int numFaces, int numEdges){

    this->numVertexes = numVertexes;
    this->numFaces = numFaces;
    this->numEdges = numEdges;
    this->numIndexes = numFaces*3;
    this->edge_data.n = 0;
    this->edge_data.a = 0;
    this->edge_data.b = 0;
    this->vertexesMalla = 0;
    this->local_mesh.v = 0;
    this->local_mesh.n = 0;
    this->local_mesh.c = 0;
    this->triangulosMalla = 0;

}
void c_malla::setCuantities( int numVertexes, int numFaces, int numEdges ){

    this->numVertexes = numVertexes;
    this->numFaces = numFaces;
    this->numEdges = numEdges;
    this->numIndexes = numFaces*3;

}

void c_malla::resetMinMax(){

    	minx=std::numeric_limits<float>::max();
	miny=std::numeric_limits<float>::max();
	minz=std::numeric_limits<float>::max();
	maxx=-1*std::numeric_limits<float>::max();
	maxy=-1*std::numeric_limits<float>::max();
	maxz=-1*std::numeric_limits<float>::max();
}

int c_malla::getNumVertices(){
    return this->numVertexes;
}

int c_malla::getNumFaces(){
    return this->numFaces;
}

int c_malla::getNumEdges(){
    return this->numEdges;
}
int c_malla::getNumIndexes(){
    return this->numIndexes;
}

float cpu_dotProduct( float3 v1, float3 v2 ){
        return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
}

float3 cpu_crossProduct( float3 v1, float3 v2){

     //Calculando N1
     return make_float3(  (v1.y * v2.z) - (v2.y * v1.z),
                        ((v1.x * v2.z) - (v2.x * v1.z)),
                         (v1.x * v2.y) - (v2.x * v1.y));
}

float cpu_magnitude( float3 v ){
        return sqrtf( powf(v.x, 2.0f) + powf(v.y, 2.0f) + powf(v.z, 2.0f));
}

//! FLOP = 6
float cpu_angleBetweenVectors( float3 v1, float3 v2 ){

    return atan2f( cpu_magnitude( cpu_crossProduct(v1, v2) ), cpu_dotProduct( v1, v2 ) );

}


//metodos CPU para hacer improve
int c_malla::cpu_2d_delaunay(int i){

    int op1 = this->triangulosMalla[this->edge_data.op[i].x]; //index vertex opuesto 1
    int op2 = this->triangulosMalla[this->edge_data.op[i].y]; //index vertex opuesto 2
    int com_a = this->triangulosMalla[this->edge_data.a[i].x];
    int com_b = this->triangulosMalla[this->edge_data.a[i].y];

    float3 aux;
    float4 pivot, TPoint;

    float4* v_data = this->local_mesh.v;

    pivot       = v_data[op1];
    TPoint      = v_data[com_a];
    aux         = make_float3(TPoint.x - pivot.x, TPoint.y - pivot.y, TPoint.z - pivot.z + 0.0*TPoint.w*pivot.w); //! + 5 flop
    TPoint      = v_data[com_b];
    float va    = cpu_angleBetweenVectors(aux, make_float3(TPoint.x - pivot.x, TPoint.y - pivot.y, TPoint.z - pivot.z + 0.0*TPoint.w*pivot.w) ); //! + 11 flop

    pivot       = v_data[op2];
    TPoint      = v_data[com_a];
    aux         = make_float3(TPoint.x - pivot.x, TPoint.y - pivot.y, TPoint.z - pivot.z + 0.0*TPoint.w*pivot.w); //! + 5 flop
    TPoint      = v_data[com_b];
    float wa    = cpu_angleBetweenVectors(aux, make_float3(TPoint.x - pivot.x, TPoint.y - pivot.y, TPoint.z - pivot.z + 0.0*TPoint.w*pivot.w)); //! + 11 flop

    return (int)(fabs(va + wa)/LAWSON_PI - EPSILON_DELAUNAY); //! + 7 flop
}

//! 62 flop
int c_malla::cpu_3d_delaunay(int i, const float limit_angle ){

    int op1 = this->triangulosMalla[this->edge_data.op[i].x]; //index vertex opuesto 1
    int op2 = this->triangulosMalla[this->edge_data.op[i].y]; //index vertex opuesto 2
    int com_a = this->triangulosMalla[this->edge_data.a[i].x];
    int com_b = this->triangulosMalla[this->edge_data.a[i].y];
    float3 aux1, aux2, n1, n2;
    float4 pivot, TPoint;
    float4* mesh_data = this->local_mesh.v;
    // get pivot, and the other two points for first triangle
    pivot       = mesh_data[op1];
    TPoint      = mesh_data[com_a];
    aux1        = make_float3(TPoint.x - pivot.x, TPoint.y - pivot.y, TPoint.z - pivot.z + 0.0*TPoint.w*pivot.w); //! + 5 flop
    TPoint      = mesh_data[com_b];
    aux2        = make_float3(TPoint.x - pivot.x, TPoint.y - pivot.y, TPoint.z - pivot.z + 0.0*TPoint.w*pivot.w); //! + 5 flop
    // compute angle
    float va    = cpu_angleBetweenVectors(aux1, aux2); //! + 6 flop
    // compute cpu_crossProduct product for 3d filter
    n1          = cpu_crossProduct(aux1, aux2); //! + 9flop

    // the same for the other triangle
    pivot       = mesh_data[op2];
    TPoint      = mesh_data[com_a];
    aux1        = make_float3(TPoint.x - pivot.x, TPoint.y - pivot.y, TPoint.z - pivot.z + 0.0*TPoint.w*pivot.w); //! + 5 flop
    TPoint      = mesh_data[com_b];
    aux2        = make_float3(TPoint.x - pivot.x, TPoint.y - pivot.y, TPoint.z - pivot.z + 0.0*TPoint.w*pivot.w); //! + 5 flop
    float wa    = cpu_angleBetweenVectors(aux1, aux2); //! + 6 flop
    n2          = cpu_crossProduct(aux2, aux1); //! + 9 flop

    float angle = fabs(cpu_angleBetweenVectors( n1, n2 ));
    return (int)(fabs(va + wa)/LAWSON_PI - EPSILON_DELAUNAY)*((int)(limit_angle/angle)); //! + 12 flop

}

void c_malla::printTriangleArray(){

    for(int i=0; i<this->numFaces; i++){
        printf("face[%i] = %i %i %i\n", i,  triangulosMalla[i*3], triangulosMalla[i*3+1], triangulosMalla[i*3+2]);
        print_edge(triangle_edges[i*3]);
        print_edge(triangle_edges[i*3+1]);
        print_edge(triangle_edges[i*3+2]);
        printf("--\n");
    }

}


int c_malla::cpu_fix_edge(int i, int rt){

    int2 n = edge_data.n[i];
    int2 a = edge_data.a[i];
    int2 b = edge_data.b[i];
    int2 op = edge_data.op[i];
    unsigned int *eab = this->triangulosMalla;
    int triangulo = 0;
    int fixed = 0;
    if( (n.x != eab[a.x] || n.y != eab[a.y]) ){

        int triangulo = rt;
        if( eab[3*triangulo+0] == n.x ){
           a.x = 3*triangulo+0;
           eab[3*triangulo+1] == n.y ? (a.y = 3*triangulo+1, op.x = 3*triangulo+2) : (a.y = 3*triangulo+2, op.x = 3*triangulo+1);
        }
        else if( eab[3*triangulo+1] == n.x ){
           a.x = 3*triangulo+1;
           eab[3*triangulo+0] == n.y ? (a.y = 3*triangulo+0, op.x = 3*triangulo+2) : (a.y = 3*triangulo+2, op.x = 3*triangulo+0);
        }
        else if( eab[3*triangulo+2] == n.x ){
           a.x = 3*triangulo+2;
           eab[3*triangulo+0] == n.y ? (a.y = 3*triangulo+0, op.x = 3*triangulo+1) : (a.y = 3*triangulo+1, op.x = 3*triangulo+0);
        }
        fixed = 1;
    }
    else if( b.x != -1 ){
        if( (n.x != eab[b.x] || n.y != eab[b.y]) ){
            int triangulo = rt;
            if( eab[3*triangulo+0] == n.x ){
               b.x = 3*triangulo+0;
               eab[3*triangulo+1] == n.y ? (b.y = 3*triangulo+1, op.y = 3*triangulo+2) : (b.y = 3*triangulo+2, op.y = 3*triangulo+1);
            }
            else if( eab[3*triangulo+1] == n.x ){
               b.x = 3*triangulo+1;
               eab[3*triangulo+0] == n.y ? (b.y = 3*triangulo+0, op.y = 3*triangulo+2) : (b.y = 3*triangulo+2, op.y = 3*triangulo+0);
            }
            else if( eab[3*triangulo+2] == n.x ){
               b.x = 3*triangulo+2;
               eab[3*triangulo+0] == n.y ? (b.y = 3*triangulo+0, op.y = 3*triangulo+1) : (b.y = 3*triangulo+1, op.y = 3*triangulo+0);
            }
            fixed = 1;
        }
    }
    edge_data.a[i] = make_int2(a.x, a.y);
    edge_data.b[i] = make_int2(b.x, b.y);
    edge_data.op[i] = make_int2(op.x, op.y);
    return fixed;

}


void c_malla::cpu_fix_add_neighbors( int i ){

    // REPARANDO LOS EDGES
    int ot = this->edge_data.a[i].x / 3;
    int ot2 = this->edge_data.b[i].x / 3;
    int e_t1_index, e_t2_index;
    for(int k=0; k<3; k++){
        int fe = triangle_edges[ot*3+k];
        if(fe != i){
            if( this->cpu_fix_edge(fe, ot2) == 1 ){
                e_t1_index = ot*3+k;
            }
            if( edge_marks[fe] == 0 ){
                this->edge_marks[fe] = 1;
                if( edge_data.b[fe].x != -1 )
                    this->edge_stack.push(fe);
            }
        }
        fe = triangle_edges[ot2*3+k];
        if(fe != i){
            // fixear y agregar edge
            if( this->cpu_fix_edge(fe, ot) == 1 ){
                e_t2_index = ot2*3+k;
            }
            if( edge_marks[fe] == 0 ){
                this->edge_marks[fe] = 1;
                if( edge_data.b[fe].x != -1 )
                    this->edge_stack.push(fe);
            }
        }
    }
    //reparando el edge_indexes_per_edge, intercambiando los dos edges que deben ser switcheados.
    int aux = this->triangle_edges[e_t1_index];
    this->triangle_edges[e_t1_index] = this->triangle_edges[e_t2_index];
    this->triangle_edges[e_t2_index] = aux;
}



void c_malla::make_stack(){

    for(int i=0; i<this->numEdges; i++){
       if(edge_data.b[i].x != -1 ){
            this->edge_stack.push(i);
            this->edge_marks[i] = 1;
        }
    }
}


void c_malla::cpu_edge_flip( int i){

    int2 n = edge_data.n[i];
    int2 a = edge_data.a[i];
    int2 b = edge_data.b[i];
    int2 op = edge_data.op[i];
    //Finalmente Intercambiar los indices necesarios
    triangulosMalla[edge_data.a[i].x] = triangulosMalla[op.y];
    triangulosMalla[edge_data.b[i].y] = triangulosMalla[op.x];
    //y por ultimo actualizar la arista modificada.
    //mantiene los mismos triangulos de referencia pero modifica sus indices a los vertices
    edge_data.a[i] = make_int2(op.x, a.x);
    edge_data.b[i] = make_int2(b.y, op.y);
    //el par b1 b2 tiene que asignar al reves cambiable2 y op1 ,
    //a diferencia del par a1 a2 que es opuesto1 y cambiable1  respectivamente
    edge_data.n[i] = make_int2(triangulosMalla[op.x], triangulosMalla[a.x]);
    edge_data.op[i] = make_int2(a.y, b.x);

}

void c_malla::print_edge(int i){

    printf("edge[%i] =  %i <----> %i    (op:  [%i]%i <-----> [%i]%i)\n", i, edge_data.n[i].x, edge_data.n[i].y,
                                        edge_data.op[i].x, triangulosMalla[edge_data.op[i].x],
                                        edge_data.op[i].y,triangulosMalla[edge_data.op[i].y]);
}


void c_malla::cpu_improve_3d(){

    int ef = 0;
    //printf("STACK COMPLETED; size = %i\n", (int)this->edge_stack.size());
    while(!this->edge_stack.empty()){
        int edge = this->edge_stack.top();  // sacar edge de la pila
        this->edge_stack.pop();             // eliminar ese elemento de la pila
        //printf("STACK = %i\n", (int)this->edge_stack.size());
        //print_edge(edge);
        this->edge_marks[edge] = 0;         // desmarcar el edge
        if (cpu_3d_delaunay(edge, 10.0f) > 0){         // consultar si necesita ser flipeado
            cpu_edge_flip(edge);            // en caso positivo, entonces flipearlo
            ef++;
            //printf("edge flipped!\n");
            cpu_fix_add_neighbors(edge);    // fixear los edges vecinos, son 4
        }
        //this->printTriangleArray();
        //printf("--------------------\n");
    }

}


int c_malla::cpu_improve(){

    int ef = 0;
    //printf("STACK COMPLETED; size = %i\n", (int)this->edge_stack.size());
    while(!this->edge_stack.empty()){
        int edge = this->edge_stack.top();  // sacar edge de la pila
        this->edge_stack.pop();             // eliminar ese elemento de la pila
        //printf("STACK = %i\n", (int)this->edge_stack.size());
        //print_edge(edge);
        this->edge_marks[edge] = 0;         // desmarcar el edge
        if (cpu_2d_delaunay(edge) > 0){         // consultar si necesita ser flipeado
            cpu_edge_flip(edge);            // en caso positivo, entonces flipearlo
            ef++;
            //printf("edge flipped!\n");
            cpu_fix_add_neighbors(edge);    // fixear los edges vecinos, son 4
        }
        //this->printTriangleArray();
        //printf("--------------------\n");
    }
    //printf("edgeflips = %i\nv=%d, e=%d, f=%i\n", ef, getNumVertices(), getNumEdges(), getNumFaces());
    return ef;
}

c_malla::~c_malla(){
    if( edge_data.n ){
        //printf("Deleting Host edge_data.n......");
        free(edge_data.n);
        edge_data.n = 0;
        //printf("OK\n");
    }
    if( edge_data.a ){
        //printf("Deleting Host edge_data.a......");
        free(edge_data.a);
        edge_data.a = 0;
        //printf("OK\n");
    }
    if( edge_data.b ){
        //printf("Deleting Host edge_data.b......");
        free(edge_data.b);
        edge_data.b = 0;
        //printf("OK\n");
    }
    if( edge_data.op ){
        //printf("Deleting Host edge_data.op......");
        free(edge_data.op);
        edge_data.op = 0;
        //printf("OK\n");
    }
    free(this->collisionsMap);
}

void c_malla::export_off(const char *filename){

	//Lineas necesarias para que scanf lea archivo en computadores seteados en otro lenguaje.
	setlocale(LC_NUMERIC, "POSIX");
	FILE *fw = fopen(filename,"w");
	fprintf(fw,"OFF\n");
	fprintf(fw,"%d %d %d\n",getNumVertices(),getNumFaces(),getNumEdges());
	for(int i=0; i<getNumVertices(); i++) {
		fprintf(fw,"%f %f %f\n",vertexesMalla[i].x,vertexesMalla[i].y,vertexesMalla[i].z);
	}
	for(int i=0; i<getNumFaces(); i++) {
		fprintf(fw,"%d %d %d %d\n", 3, triangulosMalla[i*3+0],triangulosMalla[i*3+1], triangulosMalla[i*3+2] );
	}
	fclose(fw);
	setlocale(LC_NUMERIC, "");
}


void normalize(float3 *v){
	double mag = sqrt(v->x*v->x + v->y*v->y + v->z*v->z);
	v->x /= mag;
	v->y /= mag;
	v->z /= mag;
}

double c_malla::min_angle(){
	int ignored=0;
	double avangle = 0.0;
	double mangle = LAWSON_PI;
	for(int i=0; i<getNumFaces(); i++){
		if( triangulosMalla[i*3 + 0] < 3 ||  triangulosMalla[i*3 + 1] < 3 || triangulosMalla[i*3 + 2] < 3 ){
			ignored++;
			continue;
		}
		Vertex a = vertexesMalla[triangulosMalla[i*3 + 0]];
		Vertex b = vertexesMalla[triangulosMalla[i*3 + 1]];
		Vertex c = vertexesMalla[triangulosMalla[i*3 + 2]];
		float3 v = make_float3(b.x - a.x, b.y - a.y, b.z - a.z);
		float3 w = make_float3(c.x - a.x, c.y - a.y, c.z - a.z);
		normalize(&v);
		normalize(&w);

		double val = cpu_dotProduct(v,w);
		double new_angle = acos(val);
		avangle += new_angle;
		if( new_angle < mangle ){
			mangle = new_angle;
			//printf("angle = %e rad\n", new_angle);
		}
		a = vertexesMalla[triangulosMalla[i*3 + 2]];
		b = vertexesMalla[triangulosMalla[i*3 + 0]];
		c = vertexesMalla[triangulosMalla[i*3 + 1]];
		v = make_float3(b.x - a.x, b.y - a.y, b.z - a.z);
		w = make_float3(c.x - a.x, c.y - a.y, c.z - a.z);
		normalize(&v);
		normalize(&w);

		val = cpu_dotProduct(v,w);
		new_angle = acos(val);
		avangle += new_angle;
		
		if( new_angle < mangle ){
			mangle = new_angle;
			//printf("angle = %e rad\n", new_angle);
		}
		a = vertexesMalla[triangulosMalla[i*3 + 1]];
		b = vertexesMalla[triangulosMalla[i*3 + 2]];
		c = vertexesMalla[triangulosMalla[i*3 + 0]];
		v = make_float3(b.x - a.x, b.y - a.y, b.z - a.z);
		w = make_float3(c.x - a.x, c.y - a.y, c.z - a.z);
		normalize(&v);
		normalize(&w);

		val = cpu_dotProduct(v,w);
		new_angle = acos(val);
		avangle += new_angle;
		if( new_angle < mangle ){
			mangle = new_angle;
			//printf("angle = %e rad\n", new_angle);
		}

	}
	avangle /= (double)(getNumFaces()*3);
	//printf("ignored boundary triangles = %d\n", ignored);
	return mangle;

}



void c_malla::append_result(int n, double t, int ef, const char *filename){

	FILE *fw = fopen(filename, "a");
	fprintf(fw, "%d\t\t%.5f\t\t%d\n", n, t, ef);
	fclose(fw);
}
