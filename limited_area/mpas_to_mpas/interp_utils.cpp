#include <math.h>
#include "interp_utils.h"

// Compute the great-circle distance between (lat1, lon1) and (lat2, lon2) on a
// sphere with given radius.
float sphere_distance(float lat1, float lon1, float lat2, float lon2, float radius)
{
	float arg1;

	arg1 = sqrtf( powf(sinf(0.5*(lat2-lat1)),2.0) + cosf(lat1)*cosf(lat2)*powf(sinf(0.5*(lon2-lon1)),2.0) );
	return 2.0 * radius * asinf(arg1);
}


int nearest_cell(float target_lat, float target_lon, int start_cell, int nCells,
                 int *nEdgesOnCell, int **cellsOnCell, float *latCell, float *lonCell)
{
	int iCell;
	int current_cell;
	int retval;
	float current_distance;
	float d;
	float nearest_distance;

	retval = start_cell;
	current_cell = -1;

	while (retval != current_cell) {
		current_cell = retval;
		current_distance = sphere_distance(latCell[current_cell], lonCell[current_cell], target_lat, target_lon, 1.0);
		retval = current_cell;
		nearest_distance = current_distance;
		for (int i=0; i<nEdgesOnCell[current_cell]; i++) {
			iCell = cellsOnCell[current_cell][i]-1;
			if (iCell <= nCells) {
				d = sphere_distance(latCell[iCell], lonCell[iCell], target_lat, target_lon, 1.0);
				if (d < nearest_distance) {
					retval = iCell;
					nearest_distance = d;
				}
			}
		}
	}

	return retval;
}


int nearest_vertex(float target_lat, float target_lon, int start_vertex,
                   int vertexDegree, int *nEdgesOnCell, int **verticesOnCell, int **cellsOnVertex,
                   float *latCell, float *lonCell, float *latVertex, float *lonVertex)
{
	int i, cell1, cell2, cell3, iCell, iVtx;
	int current_vertex;
	int retval;
	float cell1_dist, cell2_dist, cell3_dist;
	float current_distance, d, nearest_distance;

	retval = start_vertex;
	current_vertex = -1;
	
	while (retval != current_vertex) {
		current_vertex = retval;
		current_distance = sphere_distance(latVertex[current_vertex], lonVertex[current_vertex], target_lat, target_lon, 1.0);
		retval = current_vertex;
		nearest_distance = current_distance;
		cell1 = cellsOnVertex[current_vertex][0]-1;
		cell1_dist = sphere_distance(latCell[cell1], lonCell[cell1], target_lat, target_lon, 1.0);
		cell2 = cellsOnVertex[current_vertex][1]-1;
		cell2_dist = sphere_distance(latCell[cell2], lonCell[cell2], target_lat, target_lon, 1.0);
		if (vertexDegree == 3) {
			cell3 = cellsOnVertex[current_vertex][2]-1;
			cell3_dist = sphere_distance(latCell[cell3], lonCell[cell3], target_lat, target_lon, 1.0);
		}
		if (vertexDegree == 3) {
			if (cell1_dist < cell2_dist) {
				if (cell1_dist < cell3_dist) {
					iCell = cell1;
				}
				else {
					iCell = cell3;
				}
			}
			else {
				if (cell2_dist < cell3_dist) {
					iCell = cell2;
				}
				else {
					iCell = cell3;
				}
			}
		}
		else {
			if (cell1_dist < cell2_dist) {
				iCell = cell1;
			}
			else {
				iCell = cell2;
			}
		}
		for (int i=0; i<nEdgesOnCell[iCell]; i++) {
			iVtx = verticesOnCell[iCell][i]-1;
			d = sphere_distance(latVertex[iVtx], lonVertex[iVtx], target_lat, target_lon, 1.0);
			if (d < nearest_distance) {
				retval = iVtx;
				nearest_distance = d;
			}
		}
	}
}


// Returns the length of the great circle arc from A=(ax, ay, az) to
// B=(bx, by, bz). It is assumed that both A and B lie on the surface of the
// same sphere centered at the origin.
float mpas_arc_length(float ax, float ay, float az, float bx, float by, float bz)
{
	float r, c;
	float cx, cy, cz;

	cx = bx - ax;
	cy = by - ay;
	cz = bz - az;

	r = sqrtf(ax*ax + ay*ay + az*az);
	c = sqrtf(cx*cx + cy*cy + cz*cz);

	return 2.0 * r * asinf(c/(2.0*r));
}


inline float max(float a, float b)
{
	return (a > b) ? a : b;
}


//***********************************************************************
//
//  routine mpas_triangle_signed_area_sphere
//
//> \brief   Calculates area of a triangle on a sphere
//> \author  Matthew Hoffman
//> \date    13 January 2015
//> \details
//>  This routine calculates the area of a triangle on the surface of a sphere.
//>  Uses the spherical analog of Heron's formula.
//>  Copied from mesh generator.  A CCW winding angle is positive.
//-----------------------------------------------------------------------
float mpas_triangle_signed_area_sphere(float *a, float *b, float *c, float radius)
{
	float ab, bc, ca, semiperim, tanqe;
	float ablen[3], aclen[3], Dlen[3];

	float mpas_triangle_signed_area_sphere;
 
	ab = mpas_arc_length(a[0], a[1], a[2], b[0], b[1], b[2]) / radius;
	bc = mpas_arc_length(b[0], b[1], b[2], c[0], c[1], c[2]) / radius;
	ca = mpas_arc_length(c[0], c[1], c[2], a[0], a[1], a[2]) / radius;
	semiperim = 0.5 * (ab + bc + ca);
 
	tanqe = sqrtf(max(0.0,tanf(0.5 * semiperim) * tanf(0.5 * (semiperim - ab)) * tanf(0.5 * (semiperim - bc)) * tanf(0.5 * (semiperim - ca))));
 
	mpas_triangle_signed_area_sphere = 4.0 * radius * radius * atanf(tanqe);
 
	// computing correct signs (in similar fashion to mpas_sphere_angle)
	ablen[0] = b[0] - a[0];
	ablen[1] = b[1] - a[1];
	ablen[2] = b[2] - a[2];
 
	aclen[0] = c[0] - a[0];
	aclen[1] = c[1] - a[1];
	aclen[2] = c[2] - a[2];
 
	Dlen[0] =   (ablen[1] * aclen[2]) - (ablen[2] * aclen[1]);
	Dlen[1] = -((ablen[0] * aclen[2]) - (ablen[2] * aclen[0]));
	Dlen[2] =   (ablen[0] * aclen[1]) - (ablen[1] * aclen[0]);
 
	if ((Dlen[0]*a[0] + Dlen[1]*a[1] + Dlen[2]*a[2]) < 0.0) {
		mpas_triangle_signed_area_sphere = -1.0 * mpas_triangle_signed_area_sphere;
	}
 
	return mpas_triangle_signed_area_sphere;
}


//***********************************************************************
//
//  function mpas_wachspress_coordinates
//
//> \brief Compute the barycentric Wachspress coordinates for a polygon
//> \author  Phillip Wolfram
//> \date    01/26/2015
//> \details
//>  Computes the barycentric Wachspress coordinates for a polygon with nVertices
//>  points in R3, vertCoords for a particular pointInterp with normalized radius.
//>  Follows Gillette, A., Rand, A., Bajaj, C., 2011.
//>  Error estimates for generalized barycentric interpolation.
//>  Advances in computational mathematics 37 (3), 417–439.
//>  Optimized version of mpas_wachspress_coordinates uses optional cached B_i areas
//------------------------------------------------------------------------
void mpas_wachspress_coordinates(int nVertices, float vertCoords[][3], float *pointInterp, float *mpas_wachspress_coordinates)
{
	// computational intermediates
	float wach[nVertices];
	float wach_total; // The wachspress total weight
	int im1, i0, ip1;   // im1 = (i-1), i0 = i, ip1 = (i+1)
  
	// triangle areas to compute wachspress coordinate
	float areaA[nVertices];
	float areaB[nVertices];

	float radiusLocal;

// TODO:	radiusLocal = sqrtf(sum(vertCoords(:,1)**2))
	
	radiusLocal = 6371229.0;

	// compute areas
	for (int i=0; i<nVertices; i++) {
		// compute first area B_i
		// get vertex indices
		im1 = (nVertices + i - 1) % nVertices;
		i0  = i;
		ip1 = (nVertices + i + 1) % nVertices;
   
		// precompute B_i areas
		// always the same because B_i independent of xp,yp,zp
		// (COULD CACHE AND USE RESULT FROM ARRAY FOR FURTHER OPTIMIZATION)
		areaB[i] = mpas_triangle_signed_area_sphere(vertCoords[im1], vertCoords[i0], vertCoords[ip1], radiusLocal) / (radiusLocal*radiusLocal);
	}

	// compute areas
	for (int i=0; i<nVertices; i++) {
		// compute first area B_i
		// get vertex indices
		im1 = (nVertices + i - 1) % nVertices;
		i0  = i;
		ip1 = (nVertices + i + 1) % nVertices;
   
		// compute A_ij areas
		// must be computed each time
		areaA[i0] = mpas_triangle_signed_area_sphere(pointInterp, vertCoords[i0], vertCoords[ip1], radiusLocal) / (radiusLocal*radiusLocal);
	}

	// for each vertex compute wachpress coordinate
	for (int i=0; i<nVertices; i++) {
		wach[i] = areaB[i];
		for (int j=i+1; j<=(i + nVertices - 2); j++) {
			i0  = (nVertices + j) % nVertices;
			// accumulate products for A_ij subareas
			wach[i] = wach[i] * areaA[i0];
		}
	}

	// get summed weights for normalization
	wach_total = 0.0;
	for (int i=0; i<nVertices; i++) {
		wach_total = wach_total + wach[i];
	}

	// compute lambda
	for (int i=0; i<nVertices; i++) {
		mpas_wachspress_coordinates[i] = wach[i] / wach_total;
	}
}

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// SUBROUTINE CONVERT_LX
//
// Convert (lat,lon) to an (x, y, z) location on a sphere with specified radius.
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
void convert_lx(float *x, float *y, float *z, float radius, float lat, float lon)
{
   *z = radius * sinf(lat);
   *x = radius * cosf(lon) * cosf(lat);
   *y = radius * sinf(lon) * cosf(lat);
}

//
// Assumption: srcLevels is given in ascending order
//
void get_weights_1d(int nSrcLevels, float *srcLevels, float dstLevel, int *nSrcPts, int *srcPts, float *srcWghts)
{
	float bwght, twght;

//	std::cout << "Interpolating to " << dstLevel << std::endl;
//	for (int i=0; i<nSrcLevels; i++) {
//		std::cout << i << " " << srcLevels[i] << "\n";
//	}

	// Extrapolation below surface - use constant value
	if (dstLevel < srcLevels[0]) {
//		std::cout << "Extrapolate below the surface\n\n";
		*nSrcPts = 1;
		srcPts[0] = 0;
		srcWghts[0] = 1.0;
		return;
	}

	for (int i=1; i<nSrcLevels; i++) {
		if (dstLevel < srcLevels[i]) {
			twght = (dstLevel - srcLevels[i-1]) / (srcLevels[i] - srcLevels[i-1]);
			bwght = 1.0 - twght;
			*nSrcPts = 2;
			srcPts[0] = i-1;
			srcPts[1] = i;
			srcWghts[0] = bwght;
			srcWghts[1] = twght;
//			std::cout << "Interpolate between " << i-1 << " and " << i << " " << bwght << std::endl;
			return;
		}
	}

	// Extrapolation above top - use constant value
//	std::cout << "Extrapolate above the top\n\n";
	*nSrcPts = 1;
	srcPts[0] = nSrcLevels-1;
	srcWghts[0] = 1.0;
}
