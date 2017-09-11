#ifndef _REMAPPERCELL_H
#define _REMAPPERCELL_H

#include "RemapperBase.h"

class RemapperCell : virtual public RemapperBase {
public:
	RemapperCell();
	~RemapperCell();
	void remap(const std::type_info& t, int ndims, void *dst, void *src);
	void computeWeightsCell(int nCellsDst, int nVertLevelsSrc, int nVertLevelsDst, int vertexDegree, int *nEdgesOnCellSrc,
                                  int **verticesOnCellSrc, int **cellsOnVertexSrc,
                                  float *latCellSrc, float *lonCellSrc,
                                  float *latVertexSrc, float *lonVertexSrc,
                                  float **levelsSrc,
                                  float *latCellDst, float *lonCellDst,
                                  float **levelsDst);

private:
	//
	// Horizontal remapping fields
	//
	int nHDstPts;      // Number of horizontal destination points
	int maxHSrcPts;    // Maximum number of horizontal source points needed by any destination point
	int *nHSrcPts;     // Number of horizontal source points needed by any destination point
	int *HSrcPts;      // Source points needed for horizontal interpolation to each destination point
	int **HSrcPts2d;
	float *HSrcWghts;  // Source weights needed for horizontal interpolation to each destination point
	float **HSrcWghts2d;

	//
	// Vertical remapping fields
	//
	int nVDstPts;      // Number of vertical destination points
	int nVSrcLevels;   // Number of vertical source points
	int maxVSrcPts;    // Maximum number of vertical source points needed by any destination point
	int *nVSrcPts;     // Number of vertical source points needed by any destination point
	int **nVSrcPts2d;
	int *VSrcPts;      // Source points needed for vertical interpolation to each destination point
	int ***VSrcPts3d;
	float *VSrcWghts;  // Source weights needed for vertical interpolation to each destination point
	float ***VSrcWghts3d;

	//
	// Internal for handling remapping of 2- and 3-d fields
	//
	void remap1D(float *dst, float *src);
	void remap2D(float **dst, float **src);
	void remap3D(float ***dst, float ***src);
};
#endif
