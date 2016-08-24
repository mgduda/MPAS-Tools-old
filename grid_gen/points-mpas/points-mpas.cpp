//**************************************************
// points-mpas.cpp
//
//  Purpose:
//   
//   points-mpas.cpp is supposed to take in a triangulation defined on a sphere and a point set defined on a sphere, and create
//   a mpas grid file out of the two.
//
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license. 
// 
//  Modified:
//
//    03 December 2010
//
//  Author:
//
//    Doug Jacobsen
//
//    Modified by Michael G. Duda, 20 June 2015
//
//**************************************************


#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <iostream>
#include <fstream>
#include <tr1/unordered_set>
#include <vector>
#include <utility>
#include <math.h>
#include <assert.h>
#include <fcntl.h>
#include <unistd.h>

#include "triangulation.h"
#include "mpi.h"
#include "pnetcdf.h"

size_t min(size_t i1, size_t i2)
{
	return (i1 < i2) ? i1 : i2;
}


using namespace std;
using namespace tr1;

struct int_hasher {
	  size_t operator()(const int v) const { return v; }
};

int pt_type;
int tri_base;
double radius;
double eps;

size_t nCells, nEdges, nVertices, maxEdges, maxEdges2, TWO, vertexDegree;

int ncidp;
int ncells_dim;
int nedges_dim;
int nvertices_dim;
int maxedges_dim;
int maxedges2_dim;
int two_dim;
int vertexdegree_dim;

int xcell_var, ycell_var, zcell_var;
int latcell_var, loncell_var;
int xedge_var, yedge_var, zedge_var;
int latedge_var, lonedge_var;
int xvertex_var, yvertex_var, zvertex_var;
int latvertex_var, lonvertex_var;
int idx2cell_var, idx2edge_var, idx2vertex_var;
int cellsoncell_var;
int edgesoncell_var;
int verticesoncell_var;
int nedgesoncell_var;
int edgesonedge_var;
int cellsonedge_var;
int verticesonedge_var;
int nedgesonedge_var;
int cellsonvertex_var;
int edgesonvertex_var;
int areacell_var;
int angleedge_var;
int dcedge_var, dvedge_var;
int weightsonedge_var;
int areatriangle_var;
int kiteareasonvertex_var;
int meshdensity_var;

        
// Grid information, points, triangles, ccenters, edges
//unordered_set<pnt,pnt::hasher> edges;
unordered_set<pnt,pnt::edge_hasher> edges;
vector<pnt> edge_vec;
vector<pnt> points;
vector<pnt> ccenters;
vector<tri> triangles;

// Real connectivity arrays
vector<vector<int> > cellsOnCell, cellsOnEdge, cellsOnVertex;
vector<vector<int> > edgesOnCell, edgesOnEdge, edgesOnVertex;
vector<vector<int> > verticesOnCell, verticesOnEdge;
vector<vector<double> > weightsOnEdge;

// Unique connectivity holders
vector<unordered_set<int, int_hasher> > cellsOnCell_u;
vector<unordered_set<int, int_hasher> > cellsOnEdge_u;
vector<unordered_set<int, int_hasher> > cellsOnVertex_u;
vector<unordered_set<int, int_hasher> > edgesOnCell_u;
vector<unordered_set<int, int_hasher> > edgesOnEdge_u;
vector<unordered_set<int, int_hasher> > edgesOnVertex_u;
vector<unordered_set<int, int_hasher> > verticesOnCell_u;
vector<unordered_set<int, int_hasher> > verticesOnEdge_u;

// Grid Parameters
vector<vector<double> > kiteAreasOnVertex;
vector<double> areaCell;
vector<double> areaTriangle;
vector<double> angleEdge;
vector<double> dcEdge;
vector<double> dvEdge;

// Iterators for STL containers
vector<vector<int> >::iterator vec_int_itr;
vector<int>::iterator int_itr;
vector<vector<double> >::iterator vec_dbl_itr;
vector<double>::iterator dbl_itr;
unordered_set<int, int_hasher>::iterator u_int_itr;
vector<unordered_set<int, int_hasher> >::iterator us_itr;
unordered_set<pnt,pnt::hasher>::iterator edge_itr;
vector<pnt>::iterator point_itr;
vector<tri>::iterator tri_itr;

// Function Declarations
void readParameters();
void readPoints();
void readTriangulation();
void triangulatePoints();
void buildConnectivityArrays();
void orderConnectivityArrays();
void makeWeightsOnEdge();
int outputGridDimensions();
int outputGridAttributes();
int outputGridCoordinates();
int outputCellConnectivity();
int outputEdgeConnectivity();
int outputVertexConnectivity();
int outputCellParameters();
int outputVertexParameters();
int outputEdgeParameters();
int outputMeshDensity();
int outputVordrawArrays();
int writeGraphFile();


/********************************************************************************
 * sphere_angle
 *
 * Computes the angle between arcs AB and AC, given points A, B, and C
 * Equation numbers w.r.t. http://mathworld.wolfram.com/SphericalTrigonometry.html
 ********************************************************************************/
double sphere_angle(double ax, double ay, double az, double bx, double by, double bz, double cx, double cy, double cz)
{
	double a, b, c;		// Side lengths of spherical triangle ABC

	double ABx, ABy, ABz;	// The components of the vector AB
	double mAB;		// The magnitude of AB
	double ACx, ACy, ACz;	// The components of the vector AC
	double mAC;		// The magnitude of AC

	double Dx;		// The i-components of the cross product AB x AC
	double Dy;		// The j-components of the cross product AB x AC
	double Dz;		// The k-components of the cross product AB x AC

	double s;		// Semiperimeter of the triangle
	double sin_angle;

	a = acos(std::max(std::min(bx*cx + by*cy + bz*cz, (double)1.0), (double)-1.0));     // Eqn. (3)
	b = acos(std::max(std::min(ax*cx + ay*cy + az*cz, (double)1.0), (double)-1.0));     // Eqn. (2)
	c = acos(std::max(std::min(ax*bx + ay*by + az*bz, (double)1.0), (double)-1.0));     // Eqn. (1)

	ABx = bx - ax;
	ABy = by - ay;
	ABz = bz - az;

	ACx = cx - ax;
	ACy = cy - ay;
	ACz = cz - az;

	Dx =   (ABy * ACz) - (ABz * ACy);
	Dy = -((ABx * ACz) - (ABz * ACx));
	Dz =   (ABx * ACy) - (ABy * ACx);

	s = (double)0.5 * (a + b + c);

	sin_angle = sqrt(std::min((double)1.0 ,std::max((double)0.0 ,(sin(s-b)*sin(s-c))/(sin(b)*sin(c)))));   // Eqn. (28)

	if ((Dx*ax + Dy*ay + Dz*az) >= (double)0.0) {
		return 2.0 * asin(std::max(std::min(sin_angle, (double)1.0), (double)-1.0));
	}
	else {
		return -2.0 * asin(std::max(std::min(sin_angle, (double)1.0), (double)-1.0));
	}
}


int main(int argc, char ** argv)
{
        int ierr;

	int error_code;

        ierr = MPI_Init(&argc, &argv);

	cerr << "Reading in paramters" << endl;
	readParameters();

	cerr << " --- Points are defined on a sphere --- Radius = " << radius << endl;

	cerr << "Reading in points" << endl;
	readPoints();

	cerr << "Reading in triangles" << endl;
	readTriangulation();

	cerr << "Building connectivity arrays" << endl;
	buildConnectivityArrays();

	cerr << "Ordering connectivity arrays" << endl;
	orderConnectivityArrays();

	cerr << "Making weights on edge" << endl;
	makeWeightsOnEdge();

	cerr << "Writing grid dimensions" << endl;
	if(error_code = outputGridDimensions()){
		cerr << "Error - " << error_code << endl;
		exit(error_code);
	}
	cerr << "Writing grid attributes" << endl;
	if(error_code = outputGridAttributes()){
		cerr << "Error - " << error_code << endl;
		exit(error_code);
	}
	cerr << "Writing grid coordinates" << endl;
	if(error_code = outputGridCoordinates()){
		cerr << "Error - " << error_code << endl;
		exit(error_code);
	}
	cerr << "Writing cell connectivity" << endl;
	if(error_code = outputCellConnectivity()){
		cerr << "Error - " << error_code << endl;
		exit(error_code);
	}
	cerr << "Writing edge connectivity" << endl;
	if(error_code = outputEdgeConnectivity()){
		cerr << "Error - " << error_code << endl;
		exit(error_code);
	}
	cerr << "Writing vertex connectivity" << endl;
	if(error_code = outputVertexConnectivity()){
		cerr << "Error - " << error_code << endl;
		exit(error_code);
	}
	cerr << "Writing cell parameters" << endl;
	if(error_code = outputCellParameters()){
		cerr << "Error - " << error_code << endl;
		exit(error_code);
	}
	cerr << "Writing edge parameters" << endl;
	if(error_code = outputEdgeParameters()){
		cerr << "Error - " << error_code << endl;
		exit(error_code);
	}
	cerr << "Writing vertex parameters" << endl;
	if(error_code = outputVertexParameters()){
		cerr << "Error - " << error_code << endl;
		exit(error_code);
	}
	cerr << "Reading and writing meshDensity" << endl;
	if(error_code = outputMeshDensity()){
		cerr << "Error - " << error_code << endl;
		exit(error_code);
	}

	cerr << "Writing graph.info file" << endl;
	writeGraphFile();
	cerr << points.size() << " cells." << endl;
	cerr << edge_vec.size() << " edges." << endl;
	cerr << ccenters.size() << " vertices." << endl;

        ierr = MPI_Finalize();
	
	return 0;
}

void readParameters()
{
	ifstream params("Params");
	if(!params){
		cerr << "Params file not found. Writing default, and using default values." << endl;
		radius = 1.0;
		eps = 0.0;
		params.close();

		ofstream pout("Params");
		pout << "Is the input Cartesian or Latitude-Longitude (0 - Cartesian, 1 - Lat-lon)" << endl << "0" << endl;
		pout << "Are the triangles base zero or base one? (0 - base 0, 1 - base 1)" << endl << "0" << endl;
		pout << "What is the radius of the sphere these points are defined on?" << endl << "1.0" << endl;
		pout << "What was the convergence criteria used to make this grid?" << endl << "0.0" << endl;
		pout.close();
	} else {
		params.ignore(10000,'\n');
		params >> pt_type;
		params.ignore(10000,'\n');
		params.ignore(10000,'\n');
		params >> tri_base;
		params.ignore(10000,'\n');
		params.ignore(10000,'\n');
		params >> radius;
		params.ignore(10000,'\n');
		params.ignore(10000,'\n');
		params >> eps;	
		params.close();
	}

}

void readPoints()
{
	/******************************************************************
	 *
	 * This function reads in the point set from a file, and inserts it
	 * into a vector.
	 *
	 ******************************************************************/
	pnt p;
	double lat, lon;
	int i;
	ifstream pt_start("SaveVertices");

	i = 0;
	while(!pt_start.eof()){
		if(pt_type){
			pt_start >> lat >> lon;
			p = pntFromLatLon(lat,lon);
		} else {
			pt_start >> p;
		}
		p.idx = i;
		pt_start.ignore(10000,'\n');

		if(pt_start.good()){
			points.push_back(p);
		}
		i++;
	}
	pt_start.close();
}

void readTriangulation()
{
	/*****************************************************************
	 *
	 * This function reads in the triangulation from a file.
	 * It computes all of the circumcenters, and edges and adds them into
	 * corresponding vectors and hash tables.
	 *
	 * A hash table is used to store the edges, to ensure only one copy
	 * of each edge is kept.
	 *
	 *****************************************************************/

	pnt a, b, c, ccent, edge;
	pair<unordered_set<pnt,pnt::hasher>::iterator,bool> out_pair;
	double jv1, jv2, jv3;
	int vi1, vi2, vi3;
	int ei1, ei2, ei3;
	int i, j, junk;
	int min_vi;

	ifstream tris("SaveTriangles");

	tri t;

	i = 0;
	j = 0;
	while(!tris.eof()){
		tris >> vi1 >> vi2 >> vi3;
		ei1 = -1;
		ei2 = -1;
		ei3 = -1;

		if(tri_base == 1){
			vi1--;
			vi2--;
			vi3--;
		}

		tris.ignore(1000,'\n');

		if(tris.good()){
			a = points.at(vi1);
			b = points.at(vi2);
			c = points.at(vi3);

			vi1 = a.idx;
			vi2 = b.idx;
			vi3 = c.idx;

			if(!isCcw(a,b,c)){
				junk = vi2;
				vi2 = vi3;
				vi3 = junk;

				b = points.at(vi2);
				c = points.at(vi3);
			}

			circumcenter(a,b,c,ccent);

			ccent.normalize();

			ccent.idx = i;

			ccenters.push_back(ccent);

			edge = (a+b)/2.0;
			edge.idx = j;
			edge.isBdry = 0;
			edge.normalize();
			if(b.idx > a.idx){
				edge.vert_idx1 = a.idx;
				edge.vert_idx2 = b.idx;
			} else {
				edge.vert_idx2 = a.idx;
				edge.vert_idx1 = b.idx;
			}

			out_pair = edges.insert(edge);
			if(out_pair.second){
				edge_vec.push_back(edge);
				j++;
			}
			ei1 = (*out_pair.first).idx;

			edge = (b+c)/2.0;
			edge.idx = j;
			edge.isBdry = 0;
			edge.normalize();
			if(c.idx > b.idx){
				edge.vert_idx1 = b.idx;
				edge.vert_idx2 = c.idx;
			} else {
				edge.vert_idx2 = b.idx;
				edge.vert_idx1 = c.idx;
			}

			out_pair = edges.insert(edge);
			if(out_pair.second){
				edge_vec.push_back(edge);
				j++;
			}
			ei2 = (*out_pair.first).idx;

			edge = (c+a)/2.0;
			edge.idx = j;
			edge.isBdry = 0;
			edge.normalize();
			if(a.idx > c.idx){
				edge.vert_idx1 = c.idx;
				edge.vert_idx2 = a.idx;
			} else {
				edge.vert_idx2 = c.idx;
				edge.vert_idx1 = a.idx;
			}

			out_pair = edges.insert(edge);
			if(out_pair.second){
				edge_vec.push_back(edge);
				j++;
			}
			ei3 = (*out_pair.first).idx; 

			t = tri(vi1,vi2,vi3,i);

			t.ei1 = ei1;
			t.ei2 = ei2;
			t.ei3 = ei3;
			triangles.push_back(t);
			i++;
		}
	}

	edges.clear();
	tris.close();
}

void buildConnectivityArrays()
{
	/*************************************************************************
	 *
	 * This function takes the triangulation and point set previously read in
	 * from files, and computes the unique connectivity arrays for use in the
	 * ordering function.
	 *
	 * This is done by adding every item into a hash table, which takes care
	 * of duplicates for us. Later, we will order this data, and transfer it
	 * into a vector
	 *
	 *************************************************************************/
	pnt a, b, c;
	pnt edge1, edge2, edge3;
	double area_temp;
	double angle;
	int vi1, vi2, vi3, ei1, ei2, ei3;

	cellsOnCell.resize(points.size());
	edgesOnCell.resize(points.size());
	verticesOnCell.resize(points.size());
	cellsOnCell_u.resize(points.size());
	edgesOnCell_u.resize(points.size());
	verticesOnCell_u.resize(points.size());
	areaCell.resize(points.size());

	cellsOnEdge.resize(edge_vec.size());
	edgesOnEdge.resize(edge_vec.size());
	verticesOnEdge.resize(edge_vec.size());
	cellsOnEdge_u.resize(edge_vec.size());
	edgesOnEdge_u.resize(edge_vec.size());
	verticesOnEdge_u.resize(edge_vec.size());
	angleEdge.resize(edge_vec.size());
	dcEdge.resize(edge_vec.size());
	dvEdge.resize(edge_vec.size());

	cellsOnVertex_u.resize(ccenters.size());
	edgesOnVertex_u.resize(ccenters.size());
	cellsOnVertex.resize(ccenters.size());
	edgesOnVertex.resize(ccenters.size());
	areaTriangle.resize(ccenters.size());
	kiteAreasOnVertex.resize(ccenters.size());

	for(int i = 0; i < points.size(); i ++){
		areaCell[i] = 0.0;
	}

	for(int i = 0; i < ccenters.size(); i ++){
		areaTriangle[i] = 0.0;
		kiteAreasOnVertex[i].resize(3);
	}

	for(tri_itr = triangles.begin(); tri_itr != triangles.end(); ++tri_itr){
		vi1 = (*tri_itr).vi1;
		vi2 = (*tri_itr).vi2;
		vi3 = (*tri_itr).vi3;

		a = points.at(vi1);
		b = points.at(vi2);
		c = points.at(vi3);

		ei1 = (*tri_itr).ei1;
		ei2 = (*tri_itr).ei2;
		ei3 = (*tri_itr).ei3;

		edge1 = edge_vec.at(ei1);
		edge2 = edge_vec.at(ei2);
		edge3 = edge_vec.at(ei3);

		cellsOnCell_u[vi1].insert(vi2);
		cellsOnCell_u[vi1].insert(vi3);
		cellsOnCell_u[vi2].insert(vi3);
		cellsOnCell_u[vi2].insert(vi1);
		cellsOnCell_u[vi3].insert(vi1);
		cellsOnCell_u[vi3].insert(vi2);

		cellsOnEdge_u[ei1].insert(vi1);
		cellsOnEdge_u[ei1].insert(vi2);
		cellsOnEdge_u[ei2].insert(vi2);
		cellsOnEdge_u[ei2].insert(vi3);
		cellsOnEdge_u[ei3].insert(vi3);
		cellsOnEdge_u[ei3].insert(vi1);

		cellsOnVertex_u[(*tri_itr).idx].insert(vi1);
		cellsOnVertex_u[(*tri_itr).idx].insert(vi2);
		cellsOnVertex_u[(*tri_itr).idx].insert(vi3);

		edgesOnCell_u[vi1].insert(ei1);
		edgesOnCell_u[vi1].insert(ei3);
		edgesOnCell_u[vi2].insert(ei1);
		edgesOnCell_u[vi2].insert(ei2);
		edgesOnCell_u[vi3].insert(ei2);
		edgesOnCell_u[vi3].insert(ei3);

		edgesOnEdge_u[ei1].insert(ei2);
		edgesOnEdge_u[ei1].insert(ei3);
		edgesOnEdge_u[ei2].insert(ei1);
		edgesOnEdge_u[ei2].insert(ei3);
		edgesOnEdge_u[ei3].insert(ei1);
		edgesOnEdge_u[ei3].insert(ei2);

		edgesOnVertex_u[(*tri_itr).idx].insert(ei1);
		edgesOnVertex_u[(*tri_itr).idx].insert(ei2);
		edgesOnVertex_u[(*tri_itr).idx].insert(ei3);

		verticesOnCell_u[vi1].insert((*tri_itr).idx);
		verticesOnCell_u[vi2].insert((*tri_itr).idx);
		verticesOnCell_u[vi3].insert((*tri_itr).idx);

		verticesOnEdge_u[ei1].insert((*tri_itr).idx);
		verticesOnEdge_u[ei2].insert((*tri_itr).idx);
		verticesOnEdge_u[ei3].insert((*tri_itr).idx);

		//areaCell
		area_temp = triArea(points.at(vi1),ccenters.at((*tri_itr).idx),edge3);
		area_temp += triArea(points.at(vi1),edge1,ccenters.at((*tri_itr).idx));
		areaCell[vi1] += area_temp;
		area_temp = triArea(points.at(vi2),ccenters.at((*tri_itr).idx),edge1);
		area_temp += triArea(points.at(vi2),edge2,ccenters.at((*tri_itr).idx));
		areaCell[vi2] += area_temp;
		area_temp = triArea(points.at(vi3),ccenters.at((*tri_itr).idx),edge2);
		area_temp += triArea(points.at(vi3),edge3,ccenters.at((*tri_itr).idx));
		areaCell[vi3] += area_temp;
	}
}

void orderConnectivityArrays()
{
	/******************************************************************
	 *
	 * This function takes all of the hash tables that uniquely define the connectivity
	 * arrays, and orders them to be CCW and adjacent.
	 *
	 ******************************************************************/
	int i, j, k;

	// Order cellsOnEdge and verticesOnEdge, and make angleEdge
	j = 0;
	for(i = 0; i < edge_vec.size(); i++) {
		int cell1, cell2, vert1, vert2;
		pnt u, v, cross;
		pnt np;
		pnt edge;
		double sign;
		double ax, ay, az;
		double bx, by, bz;
		double cx, cy, cz;
		double vmag;

		angleEdge[i] = 0;
		// Ensure that u (cellsOnEdge2 - cellsOnEdge1) crossed with
		// v (verticesOnEdge2 - verticesOnEdge1) is a right handed
		// cross product
		//
		// The vectors u and v represent the perpendicular and parallel velocity
		// directions along an edge
		assert(cellsOnEdge_u[i].size() == 2);
		u_int_itr = cellsOnEdge_u[i].begin();	
		cell1 = (*u_int_itr);
		u_int_itr++;
		cell2 = (*u_int_itr);

		assert(verticesOnEdge_u[i].size() == 2);
		u_int_itr = verticesOnEdge_u[i].begin();
		vert1 = (*u_int_itr);
		u_int_itr++;
		vert2 = (*u_int_itr);

		edge = edge_vec.at(i);
		edge = gcIntersect(points.at(cell1),points.at(cell2),ccenters.at(vert1),ccenters.at(vert2));
		edge.idx = i;
		edge.isBdry = 0;
		edge_vec.at(i) = edge;

		u = points.at(cell2)-points.at(cell1);
		v = ccenters.at(vert2)-ccenters.at(vert1);

		dcEdge[i] = points.at(cell1).sphereDistance(points.at(cell2));
		dvEdge[i] = ccenters.at(vert1).sphereDistance(ccenters.at(vert2));

		cross = u.cross(v);
		sign = cross.dot(edge);
		cellsOnEdge[i].push_back(cell1);
		cellsOnEdge[i].push_back(cell2);


		if(sign < 0){
			verticesOnEdge[i].push_back(vert2);
			verticesOnEdge[i].push_back(vert1);

			v = ccenters.at(vert1)-ccenters.at(vert2);
		} else{
			verticesOnEdge[i].push_back(vert1);
			verticesOnEdge[i].push_back(vert2);
		}

		// angleEdge is either:
		// 1. The angle the positive tangential direction (v)
		//    makes with the local northward direction.
		//    or
		// 2. The angles the positive normal direction (u)
		// 	  makes with the local eastward direction.
		ax = edge.x;
		ay = edge.y;
		az = edge.z;

		if (sign >= 0.0) {
			cx = ccenters.at(vert2).x;
			cy = ccenters.at(vert2).y;
			cz = ccenters.at(vert2).z;
		}
		else {
			cx = ccenters.at(vert1).x;
			cy = ccenters.at(vert1).y;
			cz = ccenters.at(vert1).z;
		}

		vmag = sqrt((cx-ax)*(cx-ax) + (cy-ay)*(cy-ay) + (cz-az)*(cz-az));

		bx = vmag * (-cos(edge.getLon())*sin(edge.getLat())) + ax;
		by = vmag * (-sin(edge.getLon())*sin(edge.getLat())) + ay;
		bz = vmag * ( cos(edge.getLat()))                    + az;
         
		vmag = sqrt(bx*bx + by*by + bz*bz);
		bx = bx / vmag;
		by = by / vmag;
		bz = bz / vmag;

		angleEdge[i] = sphere_angle(ax, ay, az, bx, by, bz, cx, cy, cz);
	}

	//Order cellsOnVertex and edgesOnVertex, areaTriangle
	for(i = 0; i < ccenters.size(); i++) {

		/*
		 * Since all of the *OnVertex arrays should only have 3 items, it's easy to order CCW
		 *
		 * Then using the triangles, build kitesAreasOnVertex, and areaTrianlge
		 */
		pnt a, b, c;
		pnt edge1, edge2, edge3;;
		int c1, c2, c3;
		int e1, e2, e3;
		int swp_int;

		u_int_itr = cellsOnVertex_u[i].begin();
		c1 = (*u_int_itr);
		u_int_itr++;
		c2 = (*u_int_itr);
		u_int_itr++;
		c3 = (*u_int_itr);

		a = points.at(c1);
		b = points.at(c2);
		c = points.at(c3);

		if(!isCcw(a,b,c)){
			swp_int = c2;
			c2 = c3;
			c3 = swp_int;
		}

		cellsOnVertex[i].clear();
		cellsOnVertex[i].push_back(c1);
		cellsOnVertex[i].push_back(c2);
		cellsOnVertex[i].push_back(c3);

		u_int_itr = edgesOnVertex_u[i].begin();
		e1 = (*u_int_itr);
		u_int_itr++;
		e2 = (*u_int_itr);
		u_int_itr++;
		e3 = (*u_int_itr);

		edge1 = edge_vec.at(e1);
		edge2 = edge_vec.at(e2);
		edge3 = edge_vec.at(e3);

		if(!isCcw(edge1, edge2, edge3)){
			swp_int = e2;
			e2 = e3;
			e3 = swp_int;
		}

		edgesOnVertex[i].clear();
		edgesOnVertex[i].push_back(e1);
		edgesOnVertex[i].push_back(e2);
		edgesOnVertex[i].push_back(e3);

		areaTriangle[i] = triArea(points.at(c1),points.at(c2),points.at(c3));
	}

	//Order cellsOnCell, edgesOnCell and verticesOnCell
	for(i = 0; i < points.size(); i++) {
		// *
		// * Since we have the *OnEdge arrays ordered correctly, and the *OnVertex arrays ordered correctly,
		// * it should be easy to order the *OnCell arrays
		// *
		pnt cell_center;
		int cur_edge;
		int vert1, vert2;
		int found;
		size_t erased;

		// Choose a starting edge on cell.
		u_int_itr = edgesOnCell_u[i].begin();
		cur_edge = (*u_int_itr);
		cell_center = points.at(i);

		while(!edgesOnCell_u[i].empty()){
			//push this edge into edgesOnCell
			edgesOnCell[i].push_back(cur_edge);
			erased = edgesOnCell_u[i].erase(cur_edge);
			if(erased != 1){
				cerr << " Edge " << cur_edge << " not valid" << endl;
				cerr << " On cell " << cell_center << endl;
				cerr << " Available edges on cell..." << endl;
				for(u_int_itr = edgesOnCell_u[i].begin(); u_int_itr != edgesOnCell_u[i].end(); ++u_int_itr){
					cerr << (*u_int_itr) << " ";
				}
				cerr << endl;
				assert((int)erased == 1);
			}

			//Add the cell across the edge to cellsOnCell
			if(cellsOnEdge[cur_edge].at(0) == i){
				cellsOnCell[i].push_back(cellsOnEdge[cur_edge].at(1));
			} else {
				cellsOnCell[i].push_back(cellsOnEdge[cur_edge].at(0));
			}

			// Get the correct vertices on current edge, vert1 = starting vertex, vert2 = ending vertex
			// at least in the ccw order here
			vert1 = verticesOnEdge[cur_edge].at(0);
			vert2 = verticesOnEdge[cur_edge].at(1);

			if(!isCcw(cell_center, ccenters.at(vert1), ccenters.at(vert2))){
				vert1 = verticesOnEdge[cur_edge].at(1);
				vert2 = verticesOnEdge[cur_edge].at(0);
			}

			//Push the end vertex back into verticesOnCell (in the correct order).
			verticesOnCell[i].push_back(vert2);

			// Find the next edge by cycling over the edges connceted to the ending vertex. 
			found = 0;
			for(int_itr = edgesOnVertex[vert2].begin(); int_itr != edgesOnVertex[vert2].end(); ++int_itr){
				if((cellsOnEdge[(*int_itr)].at(0) == i || cellsOnEdge[(*int_itr)].at(1) == i) 
						&& ((*int_itr) != cur_edge)){
					cur_edge = (*int_itr);
					found = 1;
					break;
				}
			}

			if(found != 1){
				break;
			}
		}
	}

	//Order edgesOnEdge
	for(i = 0; i < edge_vec.size(); i++) {
		/*
		 * Edges on edge should be easily built now that all of the other connectivity arrays are built
		 */
		int cell1, cell2;
		int found;

		// Get cells connected to current edge
		cell1 = cellsOnEdge[i].at(0);
		cell2 = cellsOnEdge[i].at(1);

		// Cell 1 loops
		//Find current edge on cell 1, and add all edges from there to the end of edgesOnCell to edgesOnEdge
		//Then, add all edges from the beginning of edgesOnCell to edgesOnEdge
		//Since edgesOnCell should be CCW by now, these will all be CCW as well, and will iterate from cur_edge, around cell 1, 
		//and then back ccw around cell 2
		found = 0;
		for(int_itr = edgesOnCell.at(cell1).begin(); int_itr != edgesOnCell.at(cell1).end(); ++int_itr){
			if((*int_itr) == i){
				found = 1;
			} else if (found && (*int_itr) != i){
				edgesOnEdge[i].push_back((*int_itr));
			}
		}
		assert(found == 1);
		for(int_itr = edgesOnCell.at(cell1).begin(); int_itr != edgesOnCell.at(cell1).end(); ++int_itr){
			if((*int_itr) == i){
				found = 0;
			} else if(found && (*int_itr) != i){
				edgesOnEdge[i].push_back((*int_itr));
			}
		}
		assert(found == 0);
		//Cell 2 loops
		//Do the same thing done in the cell1 loops, just in cell 2 now.
		for(int_itr = edgesOnCell.at(cell2).begin(); int_itr != edgesOnCell.at(cell2).end(); ++int_itr){
			if((*int_itr) == i){
				found = 1;
			} else if (found && (*int_itr) != i){
				edgesOnEdge[i].push_back((*int_itr));
			}
		}
		assert(found == 1);
		for(int_itr = edgesOnCell.at(cell2).begin(); int_itr != edgesOnCell.at(cell2).end(); ++int_itr){
			if((*int_itr) == i){
				found = 0;
			} else if(found && (*int_itr) != i){
				edgesOnEdge[i].push_back((*int_itr));
			}
		}
		assert(found == 0);
	}

	cellsOnCell_u.clear();
	cellsOnEdge_u.clear();
	cellsOnVertex_u.clear();
	edgesOnCell_u.clear();
	edgesOnEdge_u.clear();
	edgesOnVertex_u.clear();
	verticesOnEdge_u.clear();
	verticesOnCell_u.clear();
}

void makeWeightsOnEdge()
{
	/************************************************************
	 *
	 * This function computes weightsOnEdge based on kiteAreasOnVertex
	 * and areaCell.
	 *
	 * The weights correspond to edgesOnEdge and allow MPAS to reconstruct the
	 * edge perpendicular (previously defined as v) velocity using the edge
	 * neighbors
	 *
	 * Weight formulation is defined in J. Thurburn, et al. JCP 2009
	 * Numerical representation of geostrophic modes on arbitrarily
	 * structured C-grids
	 *
	 * I'm not entirely sure it's correct yet
	 ************************************************************/
	int i, j, k;

	weightsOnEdge.resize(edge_vec.size());

	for(i = 0; i < edge_vec.size(); i++){
		size_t jj;
		size_t cur_edge;
		size_t prev_edge;
		double nei;
		size_t cell1, cell2;
		size_t edge1, edge2;
		size_t vert1, vert2;
		size_t neoc1, neoc2;
		double de;
		double area;
		double sum_r;

		cell1 = cellsOnEdge[i].at(0);
		cell2 = cellsOnEdge[i].at(1);

		neoc1 = edgesOnCell[cell1].size();
		neoc2 = edgesOnCell[cell2].size();

		weightsOnEdge[i].resize(edgesOnEdge[i].size());

		sum_r = 0.0;
		prev_edge = i;
		de = dcEdge.at(i);
		jj = 0;
		for(j = 0; j < neoc1-1; j++){
			cur_edge = edgesOnEdge[i].at(jj);

			//Find the vertex that is shared between prev_edge and cur_edge
			if(verticesOnEdge[prev_edge].at(0) == verticesOnEdge[cur_edge].at(0) ||
					verticesOnEdge[prev_edge].at(0) == verticesOnEdge[cur_edge].at(1)){
				vert1 = verticesOnEdge[prev_edge].at(0);
			} else if (verticesOnEdge[prev_edge].at(1) == verticesOnEdge[cur_edge].at(0) ||
					verticesOnEdge[prev_edge].at(1) == verticesOnEdge[cur_edge].at(1)){
				vert1 = verticesOnEdge[prev_edge].at(1);
			} else {
				cerr << "Edge " << prev_edge << " doesn't share a vertex with edge " << cur_edge << endl;
				exit(1);
			}
	
			//Using the vertex, cell center, prev_edge and cur_edge compute the kite area using
			//the two sub triangles, CC-PE-V and CC-V-CE
			//
			//This order is CCW
			area = triArea(points.at(cell1),edge_vec.at(prev_edge), ccenters.at(vert1));
			area += triArea(points.at(cell1), ccenters.at(vert1), edge_vec.at(cur_edge));

			for(k = 0; k < 3; k++){
				if(cellsOnVertex[vert1].at(k) == cell1){
					kiteAreasOnVertex[vert1][k] = area;
				}
			}

			//Compute running sum of area ratios for edges (using kites)
			sum_r = sum_r + area/areaCell.at(cell1);

			//Compute indicator function. -1 means inward edge normal, 1 mean outward edge normal.
			//Inward means cell center is end 0 of current edge, where as 
			//Outward mean cell center is end 1 of current edge
			if(cell1 == cellsOnEdge[cur_edge].at(0)){
				nei = 1.0;
			} else {
				nei = -1.0;
			}

			//weightsOnEdge as defined in Thuburn paper referenced above.
			//nei is indicator function, 0.5 is alpha (in equation 26)
			//sum_r is running sum of area ratios
			//dvEdge and de are the le and de terms used to scale the weights.
			weightsOnEdge[i][jj] = (0.5 - sum_r)*nei*dvEdge[cur_edge]/de;

			prev_edge = cur_edge;
			jj++;
		}

		sum_r = 0.0;
		prev_edge = i;
		for(j = 0; j < neoc2-1; j++){
			cur_edge = edgesOnEdge[i].at(jj);
			if(verticesOnEdge[prev_edge].at(0) == verticesOnEdge[cur_edge].at(0) ||
					verticesOnEdge[prev_edge].at(0) == verticesOnEdge[cur_edge].at(1)){
				vert1 = verticesOnEdge[prev_edge].at(0);
			} else if (verticesOnEdge[prev_edge].at(1) == verticesOnEdge[cur_edge].at(0) ||
					verticesOnEdge[prev_edge].at(1) == verticesOnEdge[cur_edge].at(1)){
				vert1 = verticesOnEdge[prev_edge].at(1);
			} else {
				cerr << "Edge " << prev_edge << " doesn't share a vertex with edge " << cur_edge << endl;
				cerr << "Edge " << prev_edge << " has vertices " << verticesOnEdge[prev_edge].at(0) << " " << verticesOnEdge[prev_edge].at(1) << endl;
				cerr << "Edge " << cur_edge << " has vertices " << verticesOnEdge[cur_edge].at(0) << " " << verticesOnEdge[cur_edge].at(1) << endl;
				exit(1);
			}


			area = triArea(points.at(cell2),edge_vec.at(prev_edge), ccenters.at(vert1));
			area += triArea(points.at(cell2), ccenters.at(vert1), edge_vec.at(cur_edge));

			for(k = 0; k < 3; k++){
				if(cellsOnVertex[vert1].at(k) == cell2){
					kiteAreasOnVertex[vert1][k] = area;
				}
			}

			sum_r = sum_r + area/areaCell.at(cell2);

			if(cell2 == cellsOnEdge[cur_edge].at(0)){
				nei = -1.0;
			} else {
				nei = 1.0;
			}

			weightsOnEdge[i][jj] = (0.5 - sum_r)*nei*dvEdge[cur_edge]/de;
			prev_edge = cur_edge;
			jj++;
		}
	}
}

int outputGridDimensions(void)
{
	/************************************************************************
	 *
	 * This function writes the grid dimensions to the netcdf file named
	 * outputFilename
	 *
	 * **********************************************************************/

        int ierr;

	// Return this code to the OS in case of failure.
	static const int NC_ERR = 2;
	
	// open the scvtmesh file
        ierr = ncmpi_create(MPI_COMM_WORLD, "x1.huge.nc", NC_64BIT_DATA, MPI_INFO_NULL, &ncidp);

	nCells = (size_t)points.size();

	maxEdges = 0;

	for(vec_int_itr = edgesOnCell.begin(); vec_int_itr != edgesOnCell.end(); ++vec_int_itr){
		maxEdges = std::max(maxEdges, (size_t)(*vec_int_itr).size());	
	}
	
	// check to see if the file was opened
        if(ierr != NC_NOERR) return NC_ERR;
	
	// define dimensions
	nEdges = (3 * nCells) - 6;
	nVertices = (2 * nCells) - 4;
	maxEdges2 = 2 * maxEdges;
	TWO = 2;
	vertexDegree = 3;
	
        ierr = ncmpi_def_dim(ncidp, "nCells", (MPI_Offset)nCells, &ncells_dim);
        ierr = ncmpi_def_dim(ncidp, "nEdges", (MPI_Offset)(nEdges), &nedges_dim);
        ierr = ncmpi_def_dim(ncidp, "nVertices", (MPI_Offset)(nVertices), &nvertices_dim);
        ierr = ncmpi_def_dim(ncidp, "maxEdges", (MPI_Offset)(maxEdges), &maxedges_dim);
        ierr = ncmpi_def_dim(ncidp, "maxEdges2", (MPI_Offset)(maxEdges2), &maxedges2_dim);
        ierr = ncmpi_def_dim(ncidp, "TWO", (MPI_Offset)(TWO), &two_dim);
        ierr = ncmpi_def_dim(ncidp, "vertexDegree", (MPI_Offset)(vertexDegree), &vertexdegree_dim);

	cerr << " nCells --- " << nCells << endl;
	cerr << " nEdges --- " << nEdges << " " << edge_vec.size() << endl;
	cerr << " nVertices --- " << nVertices << " " << triangles.size() << endl;
	cerr << " maxEdges --- " << maxEdges << endl;
	cerr << " maxEdges2 --- " << maxEdges2 << endl;

	return 0;
}

int outputGridAttributes(void)
{
	/************************************************************************
	 *
	 * This function writes the grid dimensions to the netcdf file named
	 * outputFilename
	 *
	 * **********************************************************************/

        int ierr;
        int ncdims[2];

	// Return this code to the OS in case of failure.
	static const int NC_ERR = 2;
	
	
	// write attributes

        ierr = ncmpi_def_var(ncidp, "latCell", NC_DOUBLE, 1, &ncells_dim, &latcell_var);
        ierr = ncmpi_def_var(ncidp, "lonCell", NC_DOUBLE, 1, &ncells_dim, &loncell_var);
        ierr = ncmpi_def_var(ncidp, "xCell", NC_DOUBLE, 1, &ncells_dim, &xcell_var);
        ierr = ncmpi_def_var(ncidp, "yCell", NC_DOUBLE, 1, &ncells_dim, &ycell_var);
        ierr = ncmpi_def_var(ncidp, "zCell", NC_DOUBLE, 1, &ncells_dim, &zcell_var);
        ierr = ncmpi_def_var(ncidp, "indexToCellID", NC_INT, 1, &ncells_dim, &idx2cell_var);
        ierr = ncmpi_def_var(ncidp, "latEdge", NC_DOUBLE, 1, &nedges_dim, &latedge_var);
        ierr = ncmpi_def_var(ncidp, "lonEdge", NC_DOUBLE, 1, &nedges_dim, &lonedge_var);
        ierr = ncmpi_def_var(ncidp, "xEdge", NC_DOUBLE, 1, &nedges_dim, &xedge_var);
        ierr = ncmpi_def_var(ncidp, "yEdge", NC_DOUBLE, 1, &nedges_dim, &yedge_var);
        ierr = ncmpi_def_var(ncidp, "zEdge", NC_DOUBLE, 1, &nedges_dim, &zedge_var);
        ierr = ncmpi_def_var(ncidp, "indexToEdgeID", NC_INT, 1, &nedges_dim, &idx2edge_var);
        ierr = ncmpi_def_var(ncidp, "latVertex", NC_DOUBLE, 1, &nvertices_dim, &latvertex_var);
        ierr = ncmpi_def_var(ncidp, "lonVertex", NC_DOUBLE, 1, &nvertices_dim, &lonvertex_var);
        ierr = ncmpi_def_var(ncidp, "xVertex", NC_DOUBLE, 1, &nvertices_dim, &xvertex_var);
        ierr = ncmpi_def_var(ncidp, "yVertex", NC_DOUBLE, 1, &nvertices_dim, &yvertex_var);
        ierr = ncmpi_def_var(ncidp, "zVertex", NC_DOUBLE, 1, &nvertices_dim, &zvertex_var);
        ierr = ncmpi_def_var(ncidp, "indexToVertexID", NC_INT, 1, &nvertices_dim, &idx2vertex_var);
        ncdims[0] = ncells_dim;
        ncdims[1] = maxedges_dim;
        ierr = ncmpi_def_var(ncidp, "cellsOnCell", NC_INT, 2, ncdims, &cellsoncell_var);
        ncdims[0] = ncells_dim;
        ncdims[1] = maxedges_dim;
        ierr = ncmpi_def_var(ncidp, "edgesOnCell", NC_INT, 2, ncdims, &edgesoncell_var);
        ncdims[0] = ncells_dim;
        ncdims[1] = maxedges_dim;
        ierr = ncmpi_def_var(ncidp, "verticesOnCell", NC_INT, 2, ncdims, &verticesoncell_var);
        ierr = ncmpi_def_var(ncidp, "nEdgesOnCell", NC_INT, 1, &ncells_dim, &nedgesoncell_var);
        ncdims[0] = nedges_dim;
        ncdims[1] = maxedges2_dim;
        ierr = ncmpi_def_var(ncidp, "edgesOnEdge", NC_INT, 2, ncdims, &edgesonedge_var);
        ncdims[0] = nedges_dim;
        ncdims[1] = two_dim;
        ierr = ncmpi_def_var(ncidp, "cellsOnEdge", NC_INT, 2, ncdims, &cellsonedge_var);
        ncdims[0] = nedges_dim;
        ncdims[1] = two_dim;
        ierr = ncmpi_def_var(ncidp, "verticesOnEdge", NC_INT, 2, ncdims, &verticesonedge_var);
        ierr = ncmpi_def_var(ncidp, "nEdgesOnEdge", NC_INT, 1, &nedges_dim, &nedgesonedge_var);
        ncdims[0] = nvertices_dim;
        ncdims[1] = vertexdegree_dim;
        ierr = ncmpi_def_var(ncidp, "cellsOnVertex", NC_INT, 2, ncdims, &cellsonvertex_var);
        ncdims[0] = nvertices_dim;
        ncdims[1] = vertexdegree_dim;
        ierr = ncmpi_def_var(ncidp, "edgesOnVertex", NC_INT, 2, ncdims, &edgesonvertex_var);
        ierr = ncmpi_def_var(ncidp, "areaCell", NC_DOUBLE, 1, &ncells_dim, &areacell_var);
        ierr = ncmpi_def_var(ncidp, "angleEdge", NC_DOUBLE, 1, &nedges_dim, &angleedge_var);
        ierr = ncmpi_def_var(ncidp, "dcEdge", NC_DOUBLE, 1, &nedges_dim, &dcedge_var);
        ierr = ncmpi_def_var(ncidp, "dvEdge", NC_DOUBLE, 1, &nedges_dim, &dvedge_var);
        ncdims[0] = nedges_dim;
        ncdims[1] = maxedges2_dim;
        ierr = ncmpi_def_var(ncidp, "weightsOnEdge", NC_DOUBLE, 2, ncdims, &weightsonedge_var);
        ierr = ncmpi_def_var(ncidp, "areaTriangle", NC_DOUBLE, 1, &nvertices_dim, &areatriangle_var);
        ncdims[0] = nvertices_dim;
        ncdims[1] = vertexdegree_dim;
        ierr = ncmpi_def_var(ncidp, "kiteAreasOnVertex", NC_DOUBLE, 2, ncdims, &kiteareasonvertex_var);
        ierr = ncmpi_def_var(ncidp, "meshDensity", NC_DOUBLE, 1, &ncells_dim, &meshdensity_var);

        ierr = ncmpi_put_att_text(ncidp, NC_GLOBAL, "on_a_sphere", (MPI_Offset)16, "YES             ");
        ierr = ncmpi_put_att_double(ncidp, NC_GLOBAL, "sphere_radius", NC_DOUBLE, (MPI_Offset)1, &radius);

        ierr = ncmpi_enddef(ncidp);
	
	return 0;
}

int outputGridCoordinates(void) 
{
	/************************************************************************
	 *
	 * This function writes the grid coordinates to the netcdf file named
	 * outputFilename
	 * This includes all cell centers, vertices, and edges.
	 * Both cartesian and lat,lon, as well as all of their indices
	 *
	 * **********************************************************************/

        MPI_Offset ncstart[1], nccount[1];
        int ierr;

	// Return this code to the OS in case of failure.
	static const int NC_ERR = 2;
	
	size_t i;
	
	double *x, *y, *z, *lat, *lon;
	int *idxTo;

	// Build and write cell coordinate arrays
	x = new double[nCells];
	y = new double[nCells];
	z = new double[nCells];
	lat = new double[nCells];
	lon = new double[nCells];
	idxTo = new int[nCells];
	i = 0;
	for(point_itr = points.begin(); point_itr != points.end(); ++point_itr){
		x[i] = (*point_itr).x;
		y[i] = (*point_itr).y;
		z[i] = (*point_itr).z;
		lat[i] = (*point_itr).getLat();
		lon[i] = (*point_itr).getLon();
		idxTo[i] = (*point_itr).idx+1;

		i++;
	}

        ierr = ncmpi_begin_indep_data(ncidp);

        ierr = ncmpi_put_var_double(ncidp, latcell_var, lat);
        ierr = ncmpi_put_var_double(ncidp, loncell_var, lon);
        ierr = ncmpi_put_var_double(ncidp, xcell_var, x);
        ierr = ncmpi_put_var_double(ncidp, ycell_var, y);
        ierr = ncmpi_put_var_double(ncidp, zcell_var, z);
        ierr = ncmpi_put_var_int(ncidp, idx2cell_var, idxTo);

	free(x);
	free(y);
	free(z);
	free(lat);
	free(lon);
	free(idxTo);
	
	//Build and write edge coordinate arrays
	x = new double[nEdges];
	y = new double[nEdges];
	z = new double[nEdges];
	lat = new double[nEdges];
	lon = new double[nEdges];
	idxTo = new int[nEdges];

	i = 0;
	for(point_itr = edge_vec.begin(); point_itr != edge_vec.end(); ++point_itr){
		x[i] = (*point_itr).x;
		y[i] = (*point_itr).y;
		z[i] = (*point_itr).z;
		lat[i] = (*point_itr).getLat();
		lon[i] = (*point_itr).getLon();
		idxTo[i] = (*point_itr).idx+1;

		i++;
	}
        ncstart[0] = 0;
        nccount[0] = min((size_t)nEdges, 22000000);
        while (ncstart[0] < (nEdges - nccount[0])) {
	        ierr = ncmpi_put_vara_double(ncidp, latedge_var, ncstart, nccount, &lat[ncstart[0]]);
        	ierr = ncmpi_put_vara_double(ncidp, lonedge_var, ncstart, nccount, &lon[ncstart[0]]);
	        ierr = ncmpi_put_vara_double(ncidp, xedge_var, ncstart, nccount, &x[ncstart[0]]);
        	ierr = ncmpi_put_vara_double(ncidp, yedge_var, ncstart, nccount, &y[ncstart[0]]);
	        ierr = ncmpi_put_vara_double(ncidp, zedge_var, ncstart, nccount, &z[ncstart[0]]);
        	ierr = ncmpi_put_vara_int(ncidp, idx2edge_var, ncstart, nccount, &idxTo[ncstart[0]]);
		ncstart[0] = ncstart[0] + nccount[0];
        }
        nccount[0] = nEdges-ncstart[0];
        if (nccount[0] > 0) {
	        ierr = ncmpi_put_vara_double(ncidp, latedge_var, ncstart, nccount, &lat[ncstart[0]]);
        	ierr = ncmpi_put_vara_double(ncidp, lonedge_var, ncstart, nccount, &lon[ncstart[0]]);
	        ierr = ncmpi_put_vara_double(ncidp, xedge_var, ncstart, nccount, &x[ncstart[0]]);
        	ierr = ncmpi_put_vara_double(ncidp, yedge_var, ncstart, nccount, &y[ncstart[0]]);
	        ierr = ncmpi_put_vara_double(ncidp, zedge_var, ncstart, nccount, &z[ncstart[0]]);
        	ierr = ncmpi_put_vara_int(ncidp, idx2edge_var, ncstart, nccount, &idxTo[ncstart[0]]);
	}


	free(x);
	free(y);
	free(z);
	free(lat);
	free(lon);
	free(idxTo);

	//Build and write vertex coordinate arrays
	x = new double[nVertices];
	y = new double[nVertices];
	z = new double[nVertices];
	lat = new double[nVertices];
	lon = new double[nVertices];
	idxTo = new int[nVertices];

	i = 0;
	for(point_itr = ccenters.begin(); point_itr != ccenters.end(); ++point_itr){
		x[i] = (*point_itr).x;
		y[i] = (*point_itr).y;
		z[i] = (*point_itr).z;
		lat[i] = (*point_itr).getLat();
		lon[i] = (*point_itr).getLon();
		idxTo[i] = (*point_itr).idx+1;

		i++;
	}
        ncstart[0] = 0;
        nccount[0] = min((size_t)nVertices, 22000000);
        while (ncstart[0] < (nVertices - nccount[0])) {
	        ierr = ncmpi_put_vara_double(ncidp, latvertex_var, ncstart, nccount, &lat[ncstart[0]]);
        	ierr = ncmpi_put_vara_double(ncidp, lonvertex_var, ncstart, nccount, &lon[ncstart[0]]);
	        ierr = ncmpi_put_vara_double(ncidp, xvertex_var, ncstart, nccount, &x[ncstart[0]]);
        	ierr = ncmpi_put_vara_double(ncidp, yvertex_var, ncstart, nccount, &y[ncstart[0]]);
	        ierr = ncmpi_put_vara_double(ncidp, zvertex_var, ncstart, nccount, &z[ncstart[0]]);
        	ierr = ncmpi_put_vara_int(ncidp, idx2vertex_var, ncstart, nccount, &idxTo[ncstart[0]]);
		ncstart[0] = ncstart[0] + nccount[0];
        }
        nccount[0] = nVertices-ncstart[0];
        if (nccount[0] > 0) {
	        ierr = ncmpi_put_vara_double(ncidp, latvertex_var, ncstart, nccount, &lat[ncstart[0]]);
        	ierr = ncmpi_put_vara_double(ncidp, lonvertex_var, ncstart, nccount, &lon[ncstart[0]]);
	        ierr = ncmpi_put_vara_double(ncidp, xvertex_var, ncstart, nccount, &x[ncstart[0]]);
        	ierr = ncmpi_put_vara_double(ncidp, yvertex_var, ncstart, nccount, &y[ncstart[0]]);
	        ierr = ncmpi_put_vara_double(ncidp, zvertex_var, ncstart, nccount, &z[ncstart[0]]);
        	ierr = ncmpi_put_vara_int(ncidp, idx2vertex_var, ncstart, nccount, &idxTo[ncstart[0]]);
	}

	free(x);
	free(y);
	free(z);
	free(lat);
	free(lon);
	free(idxTo);

	return 0;
}

int outputCellConnectivity(void) 
{
	/*****************************************************************
	 *
	 * This function writes all of the *OnCell arrays. Including
	 * cellsOnCell
	 * edgesOnCell
	 * verticesOnCell
	 * nEdgesonCell
	 *
	 * ***************************************************************/

        MPI_Offset ncstart[2], nccount[2];
        int ierr;

	// Return this code to the OS in case of failure.
	static const int NC_ERR = 2;
	
	size_t i, j;

	int *tmp_arr;
	
	// Build and write COC array
	tmp_arr = new int[nCells*maxEdges];

	for(i = 0; i < nCells; i++){
		for(j = 0; j < maxEdges; j++){
			tmp_arr[i*maxEdges + j] = 0;
		}
	}

	i = 0;
	for(vec_int_itr = cellsOnCell.begin(); vec_int_itr != cellsOnCell.end(); ++vec_int_itr){
		j = 0;
		for(int_itr = (*vec_int_itr).begin(); int_itr != (*vec_int_itr).end(); ++int_itr){
			tmp_arr[i*maxEdges + j] = (*int_itr) + 1;
			j++;
		}
		i++;
	}
cerr << "Begin writing cellsOnCell" << endl;
        ncstart[0] = 0;
        ncstart[1] = 0;
        nccount[0] = min((size_t)nCells, 22000000);
        nccount[1] = maxEdges;
        while (ncstart[0] < (nCells - nccount[0])) {
            ierr = ncmpi_put_vara_int(ncidp, cellsoncell_var, ncstart, nccount, &tmp_arr[ncstart[0]*nccount[1]]);
cerr << ierr << endl;
            ncstart[0] = ncstart[0] + nccount[0];
        }
        nccount[0] = nCells-ncstart[0];
        if (nccount[0] > 0) ierr = ncmpi_put_vara_int(ncidp, cellsoncell_var, ncstart, nccount, &tmp_arr[ncstart[0]*nccount[1]]);
cerr << "Finish writing cellsOnCell " << ierr << endl;



	// Build and write EOC array
	for(i = 0; i < nCells; i++){
		for(j = 0; j < maxEdges; j++){
			tmp_arr[i*maxEdges + j] = 0;
		}
	}

	i = 0;
	for(vec_int_itr = edgesOnCell.begin(); vec_int_itr != edgesOnCell.end(); ++vec_int_itr){
		j = 0;
		for(int_itr = (*vec_int_itr).begin(); int_itr != (*vec_int_itr).end(); ++int_itr){
			tmp_arr[i*maxEdges + j] = (*int_itr) + 1;	
			j++;
		}

		i++;
	}

cerr << "Begin writing edgesOnCell" << endl;
        ncstart[0] = 0;
        ncstart[1] = 0;
        nccount[0] = min((size_t)nCells, 22000000);
        nccount[1] = maxEdges;
        while (ncstart[0] < (nCells - nccount[0])) {
            ierr = ncmpi_put_vara_int(ncidp, edgesoncell_var, ncstart, nccount, &tmp_arr[ncstart[0]*nccount[1]]);
cerr << ierr << endl;
            ncstart[0] = ncstart[0] + nccount[0];
        }
        nccount[0] = nCells-ncstart[0];
        if (nccount[0] > 0) ierr = ncmpi_put_vara_int(ncidp, edgesoncell_var, ncstart, nccount, &tmp_arr[ncstart[0]*nccount[1]]);
cerr << "Finish writing edgesOnCell " << ierr << endl;

	// Build and write VOC array 
	for(i = 0; i < nCells; i++){
		for(j = 0; j < maxEdges; j++){
			tmp_arr[i*maxEdges + j] = 0;
		}
	}

	i = 0;
	for(vec_int_itr = verticesOnCell.begin(); vec_int_itr != verticesOnCell.end(); ++vec_int_itr){
		j = 0;
		for(int_itr = (*vec_int_itr).begin(); int_itr != (*vec_int_itr).end(); ++int_itr){
			tmp_arr[i*maxEdges + j] = (*int_itr) + 1;	
			j++;
		}
		i++;
	}

cerr << "Begin writing verticesOnCell" << endl;
        ncstart[0] = 0;
        ncstart[1] = 0;
        nccount[0] = min((size_t)nCells, 22000000);
        nccount[1] = maxEdges;
        while (ncstart[0] < (nCells - nccount[0])) {
            ierr = ncmpi_put_vara_int(ncidp, verticesoncell_var, ncstart, nccount, &tmp_arr[ncstart[0]*nccount[1]]);
cerr << ierr << endl;
            ncstart[0] = ncstart[0] + nccount[0];
        }
        nccount[0] = nCells-ncstart[0];
        if (nccount[0] > 0) ierr = ncmpi_put_vara_int(ncidp, verticesoncell_var, ncstart, nccount, &tmp_arr[ncstart[0]*nccount[1]]);
cerr << "Finish writing verticesOnCell " << ierr << endl;

	free(tmp_arr);

	//Build and write nEOC array
	tmp_arr = new int[nCells];

	i = 0;
	for(vec_int_itr = edgesOnCell.begin(); vec_int_itr != edgesOnCell.end(); ++vec_int_itr){
		tmp_arr[i] = (*vec_int_itr).size();
		i++;
	}

cerr << "Begin writing nEdgesOnCell" << endl;
//        ierr = ncmpi_put_var_int(ncidp, nedgesoncell_var, tmp_arr);
        ncstart[0] = 0;
        nccount[0] = min((size_t)nCells, 22000000);
        while (ncstart[0] < (nCells - nccount[0])) {
            ierr = ncmpi_put_vara_int(ncidp, nedgesoncell_var, ncstart, nccount, &tmp_arr[ncstart[0]]);
cerr << ierr << endl;
            ncstart[0] = ncstart[0] + nccount[0];
        }
        nccount[0] = nCells-ncstart[0];
        if (nccount[0] > 0) ierr = ncmpi_put_vara_int(ncidp, nedgesoncell_var, ncstart, nccount, &tmp_arr[ncstart[0]]);
cerr << "Finish writing nEdgesOnCell " << ierr << endl;


	verticesOnCell.clear();
	edgesOnCell.clear();

	free(tmp_arr);

	return 0;
}

int outputEdgeConnectivity(void) 
{
	/*****************************************************************
	 *
	 * This function writes all of the *OnEdge arrays. Including
	 * cellsOnEdge
	 * edgesOnEdge
	 * verticesOnEdge
	 * nEdgesOnEdge
	 *
	 * ***************************************************************/

        int ierr;

	// Return this code to the OS in case of failure.
	static const int NC_ERR = 2;
	
	size_t i, j;
 
        MPI_Offset ncstart[2], nccount[2];

	int *tmp_arr;


	// Build and write EOE array
	tmp_arr = new int[nEdges*maxEdges2];

	for(i = 0; i < nEdges; i++){
		for(j = 0; j < maxEdges2; j++){
			tmp_arr[i*maxEdges2 + j] = 0;
		}
	}


	i = 0;
	for(vec_int_itr = edgesOnEdge.begin(); vec_int_itr != edgesOnEdge.end(); ++vec_int_itr){
		j = 0;
		for(int_itr = (*vec_int_itr).begin(); int_itr != (*vec_int_itr).end(); ++int_itr){
			tmp_arr[i*maxEdges2 + j] = (*int_itr) + 1;	
			j++;
		}

		i++;
	}

cerr << "Begin writing edgesOnEdge" << endl;
        ncstart[0] = 0;
        ncstart[1] = 0;
        nccount[0] = min((size_t)nEdges, 22000000);
        nccount[1] = maxEdges2;
        while (ncstart[0] < (nEdges - nccount[0])) {
            ierr = ncmpi_put_vara_int(ncidp, edgesonedge_var, ncstart, nccount, &tmp_arr[ncstart[0]*nccount[1]]);
cerr << ierr << endl;
            ncstart[0] = ncstart[0] + nccount[0];
        }
        nccount[0] = nEdges-ncstart[0];
        if (nccount[0] > 0) ierr = ncmpi_put_vara_int(ncidp, edgesonedge_var, ncstart, nccount, &tmp_arr[ncstart[0]*nccount[1]]);
cerr << "Finish writing edgesOnEdge " << ierr << endl;


	free(tmp_arr);

	// Build and write COE array
	tmp_arr = new int[nEdges*TWO];
	for(i = 0; i < nEdges; i++){
		for(j = 0; j < TWO; j++){
			tmp_arr[i*TWO + j] = 0;
		}
	}
	i = 0;
	for(vec_int_itr = cellsOnEdge.begin(); vec_int_itr != cellsOnEdge.end(); ++vec_int_itr){
		j = 0;
		for(int_itr = (*vec_int_itr).begin(); int_itr != (*vec_int_itr).end(); ++int_itr){
			tmp_arr[i*TWO + j] = (*int_itr) + 1;	
			j++;
		}

		i++;
	}

cerr << "Begin writing cellsOnEdge" << endl;
        ncstart[0] = 0;
        ncstart[1] = 0;
        nccount[0] = min((size_t)nEdges, 22000000);
        nccount[1] = TWO;
        while (ncstart[0] < (nEdges - nccount[0])) {
            ierr = ncmpi_put_vara_int(ncidp, cellsonedge_var, ncstart, nccount, &tmp_arr[ncstart[0]*nccount[1]]);
cerr << ierr << endl;
            ncstart[0] = ncstart[0] + nccount[0];
        }
        nccount[0] = nEdges-ncstart[0];
        if (nccount[0] > 0) ierr = ncmpi_put_vara_int(ncidp, cellsonedge_var, ncstart, nccount, &tmp_arr[ncstart[0]*nccount[1]]);
cerr << "Finish writing cellsOnEdge " << ierr << endl;


	// Build VOE array
	i = 0;
	for(vec_int_itr = verticesOnEdge.begin(); vec_int_itr != verticesOnEdge.end(); ++vec_int_itr){
		j = 0;
		for(int_itr = (*vec_int_itr).begin(); int_itr != (*vec_int_itr).end(); ++int_itr){
			tmp_arr[i*TWO + j] = (*int_itr) + 1;	
			j++;
		}

		i++;
	}

cerr << "Begin writing verticesOnEdge" << endl;
        ncstart[0] = 0;
        ncstart[1] = 0;
        nccount[0] = min((size_t)nEdges, 22000000);
        nccount[1] = TWO;
        while (ncstart[0] < (nEdges - nccount[0])) {
            ierr = ncmpi_put_vara_int(ncidp, verticesonedge_var, ncstart, nccount, &tmp_arr[ncstart[0]*nccount[1]]);
cerr << ierr << endl;
            ncstart[0] = ncstart[0] + nccount[0];
        }
        nccount[0] = nEdges-ncstart[0];
        if (nccount[0] > 0) ierr = ncmpi_put_vara_int(ncidp, verticesonedge_var, ncstart, nccount, &tmp_arr[ncstart[0]*nccount[1]]);
cerr << "Finish writing verticesOnEdge " << ierr << endl;

	free(tmp_arr);

	// Build and write nEoe array
	tmp_arr = new int[nEdges];
	i = 0;
	for(vec_int_itr = edgesOnEdge.begin(); vec_int_itr != edgesOnEdge.end(); ++vec_int_itr){
		tmp_arr[i] = (*vec_int_itr).size();
		i++;
	}

cerr << "Begin writing nEdgesOnEdge" << endl;
        ncstart[0] = 0;
        nccount[0] = min((size_t)nEdges, 22000000);
        while (ncstart[0] < (nEdges - nccount[0])) {
            ierr = ncmpi_put_vara_int(ncidp, nedgesonedge_var, ncstart, nccount, &tmp_arr[ncstart[0]]);
cerr << ierr << endl;
            ncstart[0] = ncstart[0] + nccount[0];
        }
        nccount[0] = nEdges-ncstart[0];
        if (nccount[0] > 0) ierr = ncmpi_put_vara_int(ncidp, nedgesonedge_var, ncstart, nccount, &tmp_arr[ncstart[0]]);
cerr << "Finish writing nEdgesOnEdge " << ierr << endl;


// MGD can probably cut this code...
#if 0
	// Build and write bdryEdge array
	i = 0;
	for(vec_int_itr = cellsOnEdge.begin(); vec_int_itr != cellsOnEdge.end(); ++vec_int_itr){
		if((*vec_int_itr).size() != 2){
			tmp_arr[i] = 1;
		} else {
			tmp_arr[i] = 0;
		}
	}
#endif

	free(tmp_arr);

	cellsOnEdge.clear();
	edgesOnEdge.clear();
	
	return 0;
}

int outputVertexConnectivity(void) 
{
	/*****************************************************************
	 *
	 * This function writes all of the *OnVertex arrays. Including
	 * cellsOnVertex
	 * edgesOnVertex
	 *
	 * ***************************************************************/

        MPI_Offset ncstart[2], nccount[2];

        int ierr;

	// Return this code to the OS in case of failure.
	static const int NC_ERR = 2;
	
	size_t i, j;

	int *tmp_arr;

	// Build and write COV array
	tmp_arr = new int[nVertices*vertexDegree];
	
	for(i = 0; i < nVertices; i++){
		for(j = 0; j < vertexDegree; j++){
			tmp_arr[i*vertexDegree + j] = 0;
		}
	}

	i = 0;
	for(vec_int_itr = cellsOnVertex.begin(); vec_int_itr != cellsOnVertex.end(); ++vec_int_itr){
		j = 0;
		for(int_itr = (*vec_int_itr).begin(); int_itr != (*vec_int_itr).end(); ++int_itr){
			tmp_arr[i*vertexDegree + j] = (*int_itr) + 1;	
			j++;
		}
		i++;
	}


cerr << "Begin writing cellsOnVertex" << endl;
        ncstart[0] = 0;
        ncstart[1] = 0;
        nccount[0] = min((size_t)nVertices, 22000000);
        nccount[1] = vertexDegree;
        while (ncstart[0] < (nVertices - nccount[0])) {
            ierr = ncmpi_put_vara_int(ncidp, cellsonvertex_var, ncstart, nccount, &tmp_arr[ncstart[0]*nccount[1]]);
cerr << ierr << endl;
            ncstart[0] = ncstart[0] + nccount[0];
        }
        nccount[0] = nVertices-ncstart[0];
        if (nccount[0] > 0) ierr = ncmpi_put_vara_int(ncidp, cellsonvertex_var, ncstart, nccount, &tmp_arr[ncstart[0]*nccount[1]]);
cerr << "Finish writing cellsOnVertex " << ierr << endl;


	// Build and write EOV array
	for(i = 0; i < nVertices; i++){
		for(j = 0; j < vertexDegree; j++){
			tmp_arr[i*vertexDegree + j] = 0;
		}
	}
	i = 0;
	for(vec_int_itr = edgesOnVertex.begin(); vec_int_itr != edgesOnVertex.end(); ++vec_int_itr){
		j = 0;
		for(int_itr = (*vec_int_itr).begin(); int_itr != (*vec_int_itr).end(); ++int_itr){
			tmp_arr[i*vertexDegree + j] = (*int_itr) + 1;	
			j++;
		}

		i++;
	}
cerr << "Begin writing edgesOnVertex" << endl;
        ncstart[0] = 0;
        ncstart[1] = 0;
        nccount[0] = min((size_t)nVertices, 22000000);
        nccount[1] = vertexDegree;
        while (ncstart[0] < (nVertices - nccount[0])) {
            ierr = ncmpi_put_vara_int(ncidp, edgesonvertex_var, ncstart, nccount, &tmp_arr[ncstart[0]*nccount[1]]);
cerr << ierr << endl;
            ncstart[0] = ncstart[0] + nccount[0];
        }
        nccount[0] = nVertices-ncstart[0];
        if (nccount[0] > 0) ierr = ncmpi_put_vara_int(ncidp, edgesonvertex_var, ncstart, nccount, &tmp_arr[ncstart[0]*nccount[1]]);
cerr << "Finish writing edgesOnVertex " << ierr << endl;


	free(tmp_arr);

// MGD can probably cut this code...
#if 0
	// Build and write bdryVert array
	tmp_arr = new int[nVertices];
	
	i = 0;
	for(vec_int_itr = cellsOnVertex.begin(); vec_int_itr != cellsOnVertex.end(); ++vec_int_itr){
		if((*vec_int_itr).size() == vertexDegree){
			tmp_arr[i] = 0;
		} else {
			tmp_arr[i] = 1;
		}
		i++;
	}

	free(tmp_arr);
#endif

	cellsOnVertex.clear();
	edgesOnVertex.clear();

	return 0;
}

int outputCellParameters(void) 
{
	/*********************************************************
	 *
	 * This function writes all cell parameters, including
	 * 	areaCell
	 *
	 * *******************************************************/

        MPI_Offset ncstart[1], nccount[1];
        int ierr;

	// Return this code to the OS in case of failure.
	static const int NC_ERR = 2;
	
	size_t i, j;

cerr << "Begin writing areaCell" << endl;
        ncstart[0] = 0;
        nccount[0] = min((size_t)nCells, 22000000);
        while (ncstart[0] < (nCells - nccount[0])) {
            ierr = ncmpi_put_vara_double(ncidp, areacell_var, ncstart, nccount, &areaCell[ncstart[0]]);
cerr << ierr << endl;
            ncstart[0] = ncstart[0] + nccount[0];
        }
        nccount[0] = nCells-ncstart[0];
        if (nccount[0] > 0) ierr = ncmpi_put_vara_double(ncidp, areacell_var, ncstart, nccount, &areaCell[ncstart[0]]);
cerr << "Finish writing areaCell " << ierr << endl;

	areaCell.clear();

	return 0;
}

int outputEdgeParameters(void) 
{
	/*********************************************************
	 *
	 * This function writes all grid parameters, including
	 *	angleEdge
	 *	dcEdge
	 *	dvEdge
	 *	weightsOnEdge
	 *
	 * *******************************************************/

        int ierr;

	// Return this code to the OS in case of failure.
	static const int NC_ERR = 2;
	
	size_t i, j;

        MPI_Offset ncstart[2], nccount[2];

	double *tmp_arr;


	//Build and write angleEdge
cerr << "Begin writing angleEdge" << endl;
        ncstart[0] = 0;
        nccount[0] = min((size_t)nEdges, 22000000);
        while (ncstart[0] < (nEdges - nccount[0])) {
            ierr = ncmpi_put_vara_double(ncidp, angleedge_var, ncstart, nccount, &angleEdge[ncstart[0]]);
cerr << ierr << endl;
            ncstart[0] = ncstart[0] + nccount[0];
        }
        nccount[0] = nEdges-ncstart[0];
        if (nccount[0] > 0) ierr = ncmpi_put_vara_double(ncidp, angleedge_var, ncstart, nccount, &angleEdge[ncstart[0]]);
cerr << "Finish writing angleEdge " << ierr << endl;


	//Build and write dcEdge
cerr << "Begin writing dcEdge" << endl;
        ncstart[0] = 0;
        nccount[0] = min((size_t)nEdges, 22000000);
        while (ncstart[0] < (nEdges - nccount[0])) {
            ierr = ncmpi_put_vara_double(ncidp, dcedge_var, ncstart, nccount, &dcEdge[ncstart[0]]);
cerr << ierr << endl;
            ncstart[0] = ncstart[0] + nccount[0];
        }
        nccount[0] = nEdges-ncstart[0];
        if (nccount[0] > 0) ierr = ncmpi_put_vara_double(ncidp, dcedge_var, ncstart, nccount, &dcEdge[ncstart[0]]);
cerr << "Finish writing dcEdge " << ierr << endl;


	//Build and write dvEdge
cerr << "Begin writing dvEdge" << endl;
        ncstart[0] = 0;
        nccount[0] = min((size_t)nEdges, 22000000);
        while (ncstart[0] < (nEdges - nccount[0])) {
            ierr = ncmpi_put_vara_double(ncidp, dvedge_var, ncstart, nccount, &dvEdge[ncstart[0]]);
cerr << ierr << endl;
            ncstart[0] = ncstart[0] + nccount[0];
        }
        nccount[0] = nEdges-ncstart[0];
        if (nccount[0] > 0) ierr = ncmpi_put_vara_double(ncidp, dvedge_var, ncstart, nccount, &dvEdge[ncstart[0]]);
cerr << "Finish writing dvEdge " << ierr << endl;


	//Build and write weightsOnEdge
	tmp_arr = new double[nEdges*maxEdges2];

	i = 0;
	for(vec_dbl_itr = weightsOnEdge.begin(); vec_dbl_itr != weightsOnEdge.end(); ++vec_dbl_itr){
		for(j = 0; j < maxEdges2; j++){
			tmp_arr[i*maxEdges2 + j] = 0.0;
		}

		j = 0;
		for(dbl_itr = (*vec_dbl_itr).begin(); dbl_itr != (*vec_dbl_itr).end(); ++dbl_itr){
			tmp_arr[i*maxEdges2 + j] = (*dbl_itr);
			j++;
		}
		i++;
	}

cerr << "Begin writing weightsOnEdge" << endl;
        ncstart[0] = 0;
        ncstart[1] = 0;
        nccount[0] = min((size_t)nEdges, 11000000);
        nccount[1] = maxEdges2;
        while (ncstart[0] < (nEdges - nccount[0])) {
            ierr = ncmpi_put_vara_double(ncidp, weightsonedge_var, ncstart, nccount, &tmp_arr[ncstart[0]*nccount[1]]);
cerr << ierr << endl;
            ncstart[0] = ncstart[0] + nccount[0];
        }
        nccount[0] = nEdges-ncstart[0];
        if (nccount[0] > 0) ierr = ncmpi_put_vara_double(ncidp, weightsonedge_var, ncstart, nccount, &tmp_arr[ncstart[0]*nccount[1]]);
cerr << "Finish writing weightsOnEdge " << ierr << endl;

	free(tmp_arr);

	angleEdge.clear();
	dcEdge.clear();
	weightsOnEdge.clear();

	return 0;
}

int outputVertexParameters(void) 
{
	/*********************************************************
	 *
	 * This function writes all vertex parameters, including
	 * 	areaTriangle
	 * 	kiteAreasOnVertex
	 *
	 * *******************************************************/

        int ierr;

	// Return this code to the OS in case of failure.
	static const int NC_ERR = 2;
	
        size_t rets;
        int fd;
	
        MPI_Offset ncstart[2], nccount[2];

	size_t i, j;

	double *tmp_arr;

	// Build and write kiteAreasOnVertex
	tmp_arr = new double[nVertices*vertexDegree];
	i = 0;
	for(vec_dbl_itr = kiteAreasOnVertex.begin(); vec_dbl_itr != kiteAreasOnVertex.end(); ++vec_dbl_itr){
		j = 0;
		for(dbl_itr = (*vec_dbl_itr).begin(); dbl_itr != (*vec_dbl_itr).end(); ++dbl_itr){
			tmp_arr[i*vertexDegree + j] = (*dbl_itr);
			j++;
		}
		i++;
	}


cerr << "Begin writing kiteAreasOnVertex" << endl;
        ncstart[0] = 0;
        ncstart[1] = 0;
        nccount[0] = min((size_t)nVertices, 22000000);
        nccount[1] = vertexDegree;
        while (ncstart[0] < (nVertices - nccount[0])) {
            ierr = ncmpi_put_vara_double(ncidp, kiteareasonvertex_var, ncstart, nccount, &tmp_arr[ncstart[0]*nccount[1]]);
cerr << ierr << endl;
            ncstart[0] = ncstart[0] + nccount[0];
        }
        nccount[0] = nVertices-ncstart[0];
        if (nccount[0] > 0) ierr = ncmpi_put_vara_double(ncidp, kiteareasonvertex_var, ncstart, nccount, &tmp_arr[ncstart[0]*nccount[1]]);
cerr << "Finish writing kiteAreasOnVertex " << ierr << endl;


	free(tmp_arr);


	// Build and write areaTriangle
cerr << "Begin writing areaTriangle" << endl;
        ncstart[0] = 0;
        nccount[0] = min((size_t)nVertices, 22000000);
        while (ncstart[0] < (nVertices - nccount[0])) {
            ierr = ncmpi_put_vara_double(ncidp, areatriangle_var, ncstart, nccount, &areaTriangle[ncstart[0]]);
cerr << ierr << endl;
            ncstart[0] = ncstart[0] + nccount[0];
        }
        nccount[0] = nVertices-ncstart[0];
        if (nccount[0] > 0) ierr = ncmpi_put_vara_double(ncidp, areatriangle_var, ncstart, nccount, &areaTriangle[ncstart[0]]);
cerr << "Finish writing areaTriangle " << ierr << endl;

	areaTriangle.clear();
	kiteAreasOnVertex.clear();

	return 0;
}

int outputMeshDensity(void) 
{
	/***************************************************************************
	 *
	 * This function writes the meshDensity variable. Read in from the file SaveDensity
	 *
	 * *************************************************************************/

        MPI_Offset ncstart[1], nccount[1];

        int ierr;

	// Return this code to the OS in case of failure.
	static const int NC_ERR = 2;
	
	size_t i, j, k;
	int junk_int;
	double junk_dbl;

	vector<double> dbl_tmp_arr;


	//Build and write meshDensity
	dbl_tmp_arr.resize(nCells);
	ifstream celldens_in("SaveDensity");

	if (!celldens_in) {
		for(i = 0 ; i < nCells; i++) {
			dbl_tmp_arr.at(i) = 1.0;
		}
	} 
	else {
		for(i = 0; i < nCells; i++) {
			celldens_in >> dbl_tmp_arr.at(i);
			
		}
	}

	celldens_in.close();

cerr << "Begin writing meshDensity" << endl;
        ncstart[0] = 0;
        nccount[0] = min((size_t)nCells, 22000000);
        while (ncstart[0] < (nCells - nccount[0])) {
            ierr = ncmpi_put_vara_double(ncidp, meshdensity_var, ncstart, nccount, &dbl_tmp_arr.at(ncstart[0]));
cerr << ierr << endl;
            ncstart[0] = ncstart[0] + nccount[0];
        }
        nccount[0] = nCells-ncstart[0];
        if (nccount[0] > 0) ierr = ncmpi_put_vara_double(ncidp, meshdensity_var, ncstart, nccount, &dbl_tmp_arr.at(ncstart[0]));
cerr << "Finish writing meshDensity " << ierr << endl;

        ierr = ncmpi_end_indep_data(ncidp);

        ierr = ncmpi_close(ncidp);

	return 0;
}

int writeGraphFile() 
{
cerr << "Actually started to write graph file..." << endl;
	//This function writes out the graph.info file, for use with metis domain decomposition software.
	ofstream graph("graph.info");

cerr << "Going to write first line..." << endl;
	graph << points.size() << " " << edge_vec.size() << endl;

cerr << "Looping over cells..." << endl;
	for(vec_int_itr = cellsOnCell.begin(); vec_int_itr != cellsOnCell.end(); ++vec_int_itr){
		for(int_itr = (*vec_int_itr).begin(); int_itr != (*vec_int_itr).end(); ++int_itr){
			graph << (*int_itr)+1 << " ";
		}
		graph << endl;
	}

cerr << "Closing the file..." << endl;
	graph.close();
	cellsOnCell.clear();
	return 0;
}

