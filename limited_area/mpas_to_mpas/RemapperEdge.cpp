#include <iostream>
#include <typeinfo>
#include <typeindex>
#include <math.h>
#include "RemapperEdge.h"
#include "array_utils.hpp"
#include "interp_utils.h"

RemapperEdge::RemapperEdge()
{
	nHDstPts = 0;
	maxHSrcPts = 0;
	nHSrcPts = NULL;
	HSrcPts = NULL;
	HSrcPts2d = NULL;
	HSrcWghts = NULL;
	HSrcWghts2d = NULL;

	nVDstPts = 0;
	nVSrcLevels = 0;
	maxVSrcPts = 0;
	nVSrcPts = NULL;
	nVSrcPts2d = NULL;
	VSrcPts = NULL;
	VSrcPts3d = NULL;
	VSrcWghts = NULL;
	VSrcWghts3d = NULL;
}

RemapperEdge::~RemapperEdge()
{
	if (nHSrcPts != NULL) delete[] nHSrcPts;
	if (HSrcPts != NULL) delete[] HSrcPts;
	if (HSrcWghts != NULL) delete[] HSrcWghts;

	if (HSrcPts2d != NULL) {
		deallocate_2d<int>(HSrcPts2d);
	}
	if (HSrcWghts2d != NULL) {
		deallocate_2d<float>(HSrcWghts2d);
	}

	if (nVSrcPts != NULL) delete[] nVSrcPts;
	if (VSrcPts != NULL) delete[] VSrcPts;
	if (VSrcWghts != NULL) delete[] VSrcWghts;

	if (nVSrcPts2d != NULL) {
		deallocate_2d<int>(nVSrcPts2d);
	}
	if (VSrcPts3d != NULL) {
		deallocate_3d<int>(VSrcPts3d);
	}
	if (VSrcWghts3d != NULL) {
		deallocate_3d<float>(VSrcWghts3d);
	}
}

void RemapperEdge::computeWeightsEdge(int nCellsSrc, int nEdgesDst, int nVertLevelsSrc, int nVertLevelsDst,
                                          int *nEdgesOnCellSrc, int **cellsOnCellSrc, int **edgesOnCellSrc,
                                          float *latCellSrc, float *lonCellSrc,
                                          float *latEdgeSrc, float *lonEdgeSrc,
                                          float **levelsSrc,
                                          float *latCellDst, float *lonCellDst,
                                          float *latEdgeDst, float *lonEdgeDst,
                                          float **levelsDst)
{
	int j;
	const int maxEdges = 7;      // TODO: pass this from source mesh
	float tempLevels[nVertLevelsSrc];
	float vertCoords[maxEdges][3];
	float pointInterp[3];

	nHDstPts = nEdgesDst;
	maxHSrcPts = maxEdges;
	nHSrcPts = new int[nHDstPts];   // set to all 1 for now...
	HSrcPts = new int[nHDstPts * maxHSrcPts];
	HSrcPts2d = allocate_2d<int>(nHDstPts, maxHSrcPts, HSrcPts);
	HSrcWghts = new float[nHDstPts * maxHSrcPts];
	HSrcWghts2d = allocate_2d<float>(nHDstPts, maxHSrcPts, HSrcWghts);

	j = 0;
#pragma omp parallel for firstprivate(j) private(pointInterp, vertCoords)
	for (int i=0; i<nHDstPts; i++) {
		j = nearest_cell(latEdgeDst[i], lonEdgeDst[i], j, nCellsSrc,
                                 nEdgesOnCellSrc, cellsOnCellSrc, latCellSrc, lonCellSrc);
		nHSrcPts[i] = nEdgesOnCellSrc[j];

		convert_lx(&pointInterp[0], &pointInterp[1], &pointInterp[2], 6371229.0, latEdgeDst[i], lonEdgeDst[i]);
		for (int k=0; k<nHSrcPts[i]; k++) {
			HSrcPts2d[i][k] = edgesOnCellSrc[j][k] - 1;
			convert_lx(&vertCoords[k][0], &vertCoords[k][1], &vertCoords[k][2], 6371229.0, latEdgeSrc[HSrcPts2d[i][k]], lonEdgeSrc[HSrcPts2d[i][k]]);
		}

		mpas_wachspress_coordinates(nHSrcPts[i], vertCoords, pointInterp, HSrcWghts2d[i]);
	}

	if (nVertLevelsSrc > 0 && nVertLevelsDst > 0) {
		nVDstPts = nVertLevelsDst;
		nVSrcLevels = nVertLevelsSrc;
		maxVSrcPts = 2;
		nVSrcPts = new int[nHDstPts * nVDstPts];   // set to all 1 for now...
		nVSrcPts2d = allocate_2d<int>(nHDstPts, nVDstPts, nVSrcPts);
		VSrcPts = new int[nHDstPts * nVDstPts * maxVSrcPts];
		VSrcPts3d = allocate_3d<int>(nHDstPts, nVDstPts, maxVSrcPts, VSrcPts);
		VSrcWghts = new float[nHDstPts * nVDstPts * maxVSrcPts];
		VSrcWghts3d = allocate_3d<float>(nHDstPts, nVDstPts, maxVSrcPts, VSrcWghts);

#pragma omp parallel for private(tempLevels)
		for (int i=0; i<nHDstPts; i++) {
			// Horizontally interpolate column of levelsSrc values
			for (int k=0; k<nVertLevelsSrc; k++) {
				tempLevels[k] = 0;
			}
// TODO: How to handle the fact that zgrid is on cells? maybe just pre-compute zgridEdge?
			for (int j=0; j<nHSrcPts[i]; j++) {
				for (int k=0; k<nVertLevelsSrc; k++) {
					tempLevels[k] += (HSrcWghts2d[i][j] * levelsSrc[HSrcPts2d[i][j]][k]);
				}
			}

			// For each vertical destination point, determine weights from tempLevels
			for (int k=0; k<nVDstPts; k++) {
				get_weights_1d(nVertLevelsSrc, tempLevels, levelsDst[i][k], &nVSrcPts2d[i][k], VSrcPts3d[i][k], VSrcWghts3d[i][k]);
			}
		}
	}
}

void RemapperEdge::remap(const std::type_info& t, int ndims, void *dst, void *src)
{
	if (std::type_index(t) == typeid(float)) {
		if (ndims == 1) {
			float *dstf = (float *)dst;
			float *srcf = (float *)src;
			remap1D(dstf, srcf);
		}
		else if (ndims == 2) {
			float **dstf = (float **)dst;
			float **srcf = (float **)src;
			remap2D(dstf, srcf);
		}
		else if (ndims == 3) {
			float ***dstf = (float ***)dst;
			float ***srcf = (float ***)src;
			remap3D(dstf, srcf);
		}
	}
	else {
		throw "RemapperEdge can only handle 2-d or 3-d float fields";
	}
}

void RemapperEdge::remap1D(float *dst, float *src)
{
	std::cerr << "Remapping 1d field\n";

	for (int i=0; i<nHDstPts; i++) {
		dst[i] = 0;
		for (int j=0; j<nHSrcPts[i]; j++) {
			dst[i] += (HSrcWghts2d[i][j] * src[HSrcPts2d[i][j]]);
		}
	}
}

void RemapperEdge::remap2D(float **dst, float **src)
{
	std::cerr << "Remapping 2d field\n";

	// TODO: Right now, the time dimension is the first dimension
	for (int i=0; i<nHDstPts; i++) {
		dst[0][i] = 0;
		for (int j=0; j<nHSrcPts[i]; j++) {
			dst[0][i] += (HSrcWghts2d[i][j] * src[0][HSrcPts2d[i][j]]);
		}
	}
}

void RemapperEdge::remap3D(float ***dst, float ***src)
{
	std::cerr << "Remapping 3d field\n";

	float tempLevels[nVSrcLevels];

	// TODO: Right now, the time dimension is the first dimension
#pragma omp parallel for private(tempLevels)
	for (int i=0; i<nHDstPts; i++) {
		// Horizontally interpolate column of levelsSrc values
		for (int k=0; k<nVSrcLevels; k++) {
			tempLevels[k] = 0;
		}
		for (int j=0; j<nHSrcPts[i]; j++) {
			for (int k=0; k<nVSrcLevels; k++) {
				tempLevels[k] += (HSrcWghts2d[i][j] * src[0][HSrcPts2d[i][j]][k]);
			}
		}
	
		// For each vertical destination point, interpolate
		for (int k=0; k<nVDstPts; k++) {
			dst[0][i][k] = 0;
			for (int j=0; j<nVSrcPts2d[i][k]; j++) {
				dst[0][i][k] += VSrcWghts3d[i][k][j] * tempLevels[VSrcPts3d[i][k][j]];
			}
		}
	}
}
