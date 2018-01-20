#ifndef _REMAPPEREDGE_H
#define _REMAPPEREDGE_H

#include "RemapperBase.h"

class RemapperEdge : virtual public RemapperBase {
public:
	RemapperEdge();
	~RemapperEdge();
	void remap(const std::type_info& t, int ndims, void *dst, void *src);
	void computeWeightsEdge(int nCellsSrc, int nEdgesDst, int nVertLevelsSrc, int nVertLevelsDst,
                                  int *nEdgesOnCellSrc, int **cellsOnCellSrc, int **edgesOnCellSrc,
                                  float *latCellSrc, float *lonCellSrc,
                                  float *latEdgeSrc, float *lonEdgeSrc,
                                  float **levelsSrc,
                                  float *latCellDst, float *lonCellDst,
                                  float *latEdgeDst, float *lonEdgeDst,
                                  float **levelsDst,
                                  int *maskDst);

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
	unsigned char *HDstMask;     // Mask for horizontal destination points

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
