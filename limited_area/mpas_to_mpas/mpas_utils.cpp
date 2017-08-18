#include <iostream>
#include <math.h>
#include <omp.h>


void uncouple(int nVertLevels, int nCells, float **field, float **metric)
{
#pragma omp parallel for
	for (int i=0; i<nCells; i++) {
		for (int k=0; k<nVertLevels; k++) {
			field[i][k] *= metric[i][k];
		}
	}
}


void couple(int nVertLevels, int nCells, float **field, float **metric)
{
#pragma omp parallel for
	for (int i=0; i<nCells; i++) {
		for (int k=0; k<nVertLevels; k++) {
			field[i][k] /= metric[i][k];
		}
	}
}


void reconstruct_v(int nEdges, int nVertLevels, int *nEdgesOnEdge, int **edgesOnEdge, float **weightsOnEdge, float **u, float **v)
{
	int i, k, j, jj;
	std::cerr << "Reconstructing v field for " << nEdges << " edges\n";

#pragma omp parallel for private(k,j,jj)
	for (i=0; i<nEdges; i++) {
		for (k=0; k<nVertLevels; k++) {
			v[i][k] = 0.0;
		}
		for (j=0; j<nEdgesOnEdge[i]; j++) {
			jj = edgesOnEdge[i][j] - 1;   // edgesOnEdge field is 1-based
			for (k=0; k<nVertLevels; k++) {
				v[i][k] += u[jj][k] * weightsOnEdge[i][j];
			}
		}
	}
}

void rotate_winds(int nEdges, int nVertLevels, float *angleEdge, float **u, float **v, int toEarth)
{
	int i, k;
	float uearth, vearth;
	std::cerr << "Rotating {u,v} field for " << nEdges << " edges\n";

	//
	//  Rotate to earth-relative winds
	//
	if (toEarth != 0) {
#pragma omp parallel for private(k,uearth,vearth)
		for (i=0; i<nEdges; i++) {
			for (k=0; k<nVertLevels; k++) {
				uearth = u[i][k] * cosf(angleEdge[i]) - v[i][k] * sinf(angleEdge[i]);
				vearth = v[i][k] * cosf(angleEdge[i]) + u[i][k] * sinf(angleEdge[i]);
				u[i][k] = uearth;
				v[i][k] = vearth;
			}
		}
	}
	//
	//  Rotate to normal/tangential winds
	//
	else {
#pragma omp parallel for private(k,uearth,vearth)
		for (i=0; i<nEdges; i++) {
			for (k=0; k<nVertLevels; k++) {
				uearth = u[i][k] * cosf(angleEdge[i]) + v[i][k] * sinf(angleEdge[i]);
				vearth = v[i][k] * cosf(angleEdge[i]) - u[i][k] * sinf(angleEdge[i]);
				u[i][k] = uearth;
				v[i][k] = vearth;
			}
		}
	}
}

void avg_to_midpoint(int nCells, int nLevels, float **levels, float **layers)
{
	int i, k;

#pragma omp parallel for private(k)
	for (i=0; i<nCells; i++) {
		for (k=0; k<nLevels-1; k++) {
			layers[i][k] = 0.5 * (levels[i][k] + levels[i][k+1]);
		}
	}
}

void avg_cell_to_edge(int nEdges, int nLevels, int **cellsOnEdge, float **cellfield, float **edgefield)
{
	int i, j, k;

#pragma omp parallel for private(k,j)
	for (i=0; i<nEdges; i++) {
		if (cellsOnEdge[i][0] > 0 && cellsOnEdge[i][1] > 0) {
			for (k=0; k<nLevels; k++) {
				edgefield[i][k] = 0.5 * (cellfield[cellsOnEdge[i][0]-1][k] + cellfield[cellsOnEdge[i][1]-1][k] );
			}
		}
		// TODO: We should probably just get zedge from the limited-area mesh
		else {
			if (cellsOnEdge[i][0] > cellsOnEdge[i][1]) {
				j = cellsOnEdge[i][0];
			}
			else {
				j = cellsOnEdge[i][1];
			}
			for (k=0; k<nLevels; k++) {
				edgefield[i][k] = cellfield[j-1][k];
			}
		}
	}
}
