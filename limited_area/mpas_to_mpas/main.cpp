#include <iostream>
#include <stdio.h>
#include <string.h>
#include "netcdf.h"
#include "NCField.hpp"
#include "RemapperCell.h"
#include "RemapperEdge.h"
#include "mpas_utils.h"

void start_timer(int n);
void stop_timer(int n, int *secs, int *n_secs);


int main(int argc, char **argv)
{
	int ncid;
	int stat;

	const int NUM_SCALARS = 6;
	char qxNames[NUM_SCALARS][3] = {"qv", "qc", "qr", "qi", "qs", "qg"};

	NCField<float> *latCellDst;
	NCField<float> *lonCellDst;
	NCField<float> *latEdgeDst;
	NCField<float> *lonEdgeDst;
	NCField<float> *angleEdgeDst;
	NCField<int> *cellsOnEdgeDst;
	NCField<float> *zgridDst;
	NCField<float> *zmidDst;
	NCField<float> *zedgeDst;
	NCField<float> *uDst;
	NCField<float> *vDst;
	NCField<float> *theta_mDst;
	NCField<float> *rho_zzDst;
	NCField<float> *zzDst;
	NCField<float> *zzEdgeDst;
	NCField<float> *rho_edgeDst;
	NCField<float> *wDst;
	NCField<float> *qxDst[NUM_SCALARS];
	NCField<char> *xtime;
	float ***uDstArr;
	float ***vDstArr;
	float **zgridDstArr;
	float **zmidDstArr;
	float **zzDstArr;
	float **zzEdgeDstArr;
	float ***rho_zzDstArr;
	float ***rho_edgeDstArr;
	float *angleEdgeDstArr;

	char qxLbcName[64];

	const char *globalMeshFile;
	const char *regionalMeshFile;
	const char *globalFieldFile;
	char regionalFieldFile[64];
	char **xtimeArr;
	char date[14];

	NCField<float> *latCellSrc;
	NCField<float> *lonCellSrc;
	NCField<float> *latEdgeSrc;
	NCField<float> *lonEdgeSrc;
	NCField<float> *latVertexSrc;
	NCField<float> *lonVertexSrc;
	NCField<int> *nEdgesOnCellSrc;
	NCField<int> *nEdgesOnEdgeSrc;
	NCField<int> *cellsOnCellSrc;
	NCField<int> *verticesOnCellSrc;
	NCField<int> *cellsOnVertexSrc;
	NCField<int> *cellsOnEdgeSrc;
	NCField<int> *edgesOnCellSrc;
	NCField<int> *edgesOnEdgeSrc;
	NCField<float> *weightsOnEdgeSrc;
	NCField<float> *angleEdgeSrc;
	NCField<float> *zgridSrc;
	NCField<float> *zmidSrc;
	NCField<float> *zedgeSrc;
	NCField<float> *uSrc;
	NCField<float> *vSrc;
	NCField<float> *theta_mSrc;
	NCField<float> *rho_zzSrc;
	NCField<float> *zzSrc;
	NCField<float> *wSrc;
	NCField<float> *qxSrc[NUM_SCALARS];
	float ***uSrcArr;
	float ***vSrcArr;
	float **zgridSrcArr;
	float **zmidSrcArr;
	float **zzSrcArr;
	float ***rho_zzSrcArr;
	int *nEdgesOnEdgeSrcArr;
	int **edgesOnEdgeSrcArr;
	float **weightsOnEdgeSrcArr;
	float *angleEdgeSrcArr;
	RemapperCell *cellLayerMap;
	RemapperCell *cellLevelMap;
	RemapperCell *cellToEdgeMap;
	RemapperEdge *edgeMap;
	int secs, nsecs;
	int itime;
	int argv_idx;
	int use_reconstruct_winds;


	if (argc < 4) {
		std::cerr << "\nUsage: " << argv[0] << " [--use-reconstruct-winds] <global_IC_file> <regional_IC_file> <global_fields_file>+\n\n";
		return 1;
	}

	if (strcmp(argv[1], "--use-reconstruct-winds") == 0) {
		std::cout << "Using reconstructed winds at cell centers...\n";
		argv_idx = 2;
		use_reconstruct_winds = 1;
	}
	else {
		//
		// Try to catch other uses options besides --use-reconstruct-winds
		//
		if (argv[1][0] == '-') {
			std::cerr << "Unrecognized option: " << argv[1] << std::endl;
			std::cerr << "\nSupported options are: --use-reconstruct-winds\n";
			return 1;
		}

		std::cout << "Using normal component of winds at edges...\n";
		argv_idx = 1;
		use_reconstruct_winds = 0;
	}
	globalMeshFile = argv[argv_idx++];
	regionalMeshFile = argv[argv_idx++];


	//
	// Read mesh description fields from regional IC file
	//
	start_timer(0);
	try {
		latCellDst = new NCField<float>(regionalMeshFile, "latCell");
	}
	catch (int e) {
		std::cerr << "Error reading latCell field from " << globalMeshFile << std::endl;
		return 1;
	}
	lonCellDst = new NCField<float>(regionalMeshFile, "lonCell");
	latEdgeDst = new NCField<float>(regionalMeshFile, "latEdge");
	lonEdgeDst = new NCField<float>(regionalMeshFile, "lonEdge");
	angleEdgeDst = new NCField<float>(regionalMeshFile, "angleEdge");
	cellsOnEdgeDst = new NCField<int>(regionalMeshFile, "cellsOnEdge");
	zgridDst = new NCField<float>(regionalMeshFile, "zgrid");
	zzDst = new NCField<float>(regionalMeshFile, "zz");
	stop_timer(0, &secs, &nsecs);
	printf("Time to read mesh fields from %s : %i.%9.9i\n", regionalMeshFile, secs, nsecs);


	//
	// Read mesh description fields from global IC file
	//
	start_timer(0);
	try {
		latCellSrc = new NCField<float>(globalMeshFile, "latCell");
	}
	catch (int e) {
		std::cerr << "Error reading latCell field from " << globalMeshFile << std::endl;
		return 1;
	}
	lonCellSrc = new NCField<float>(globalMeshFile, "lonCell");
	latEdgeSrc = new NCField<float>(globalMeshFile, "latEdge");
	lonEdgeSrc = new NCField<float>(globalMeshFile, "lonEdge");
	latVertexSrc = new NCField<float>(globalMeshFile, "latVertex");
	lonVertexSrc = new NCField<float>(globalMeshFile, "lonVertex");
	nEdgesOnCellSrc = new NCField<int>(globalMeshFile, "nEdgesOnCell");
	nEdgesOnEdgeSrc = new NCField<int>(globalMeshFile, "nEdgesOnEdge");
	weightsOnEdgeSrc = new NCField<float>(globalMeshFile, "weightsOnEdge");
	edgesOnEdgeSrc = new NCField<int>(globalMeshFile, "edgesOnEdge");
	angleEdgeSrc = new NCField<float>(globalMeshFile, "angleEdge");
	cellsOnCellSrc = new NCField<int>(globalMeshFile, "cellsOnCell");
	verticesOnCellSrc = new NCField<int>(globalMeshFile, "verticesOnCell");
	cellsOnVertexSrc = new NCField<int>(globalMeshFile, "cellsOnVertex");
	cellsOnEdgeSrc = new NCField<int>(globalMeshFile, "cellsOnEdge");
	edgesOnCellSrc = new NCField<int>(globalMeshFile, "edgesOnCell");
	zgridSrc = new NCField<float>(globalMeshFile, "zgrid");
	zzSrc = new NCField<float>(globalMeshFile, "zz");
	stop_timer(0, &secs, &nsecs);
	printf("Time to read mesh fields from %s : %i.%9.9i\n", globalMeshFile, secs, nsecs);


	//
	// Create fields to hold zgrid averaged to layer midpoints and then averaged to edges
	// for both the regional and the global meshes
	//
	zmidSrc = new NCField<float>("zmid", 2, "nCells", zgridSrc->dimSize("nCells"), "nVertLevels", zgridSrc->dimSize("nVertLevelsP1")-1);
	zedgeSrc = new NCField<float>("zedge", 2, "nEdges", latEdgeSrc->dimSize("nEdges"), "nVertLevels", zgridSrc->dimSize("nVertLevelsP1")-1);
	zmidDst = new NCField<float>("zmid", 2, "nCells", zgridDst->dimSize("nCells"), "nVertLevels", zgridDst->dimSize("nVertLevelsP1")-1);
	zedgeDst = new NCField<float>("zedge", 2, "nEdges", latEdgeDst->dimSize("nEdges"), "nVertLevels", zgridDst->dimSize("nVertLevelsP1")-1);
	zzEdgeDst = new NCField<float>("zzedge", 2, "nEdges", latEdgeDst->dimSize("nEdges"), "nVertLevels", zedgeDst->dimSize("nVertLevels"));


	//
	// Average global grid zgrid field to layer midpoints and then to edges
	//
	start_timer(0);
	zmidSrcArr = zmidSrc->ptr2D();
	zgridSrcArr = zgridSrc->ptr2D();
	avg_to_midpoint((int)zmidSrc->dimSize("nCells"), (int)zmidSrc->dimSize("nVertLevels")+1, zgridSrcArr, zmidSrcArr);
	avg_cell_to_edge((int)latEdgeSrc->dimSize("nEdges"), (int)zmidSrc->dimSize("nVertLevels"), cellsOnEdgeSrc->ptr2D(), zmidSrc->ptr2D(), zedgeSrc->ptr2D());
	stop_timer(0, &secs, &nsecs);
	printf("Time to average zgridSrc to edges: %i.%9.9i\n", secs, nsecs);


	//
	// Average regional grid zgrid field to layer midpoints and then to edges
	//
	start_timer(0);
	zmidDstArr = zmidDst->ptr2D();
	zgridDstArr = zgridDst->ptr2D();
	avg_to_midpoint((int)zmidDst->dimSize("nCells"), (int)zmidDst->dimSize("nVertLevels")+1, zgridDstArr, zmidDstArr);
	avg_cell_to_edge((int)latEdgeDst->dimSize("nEdges"), (int)zmidDst->dimSize("nVertLevels"), cellsOnEdgeDst->ptr2D(), zmidDst->ptr2D(), zedgeDst->ptr2D());
	stop_timer(0, &secs, &nsecs);
	printf("Time to average zgridDst to edges: %i.%9.9i\n", secs, nsecs);


	//
	// Compute weights for remapping cell-based fields on levels (only used for w right now)
	//
	start_timer(0);
	cellLevelMap = new RemapperCell();
	cellLevelMap->computeWeightsCell(latCellDst->dimSize("nCells"), zgridSrc->dimSize("nVertLevelsP1"), zgridDst->dimSize("nVertLevelsP1"), 3,
                                      nEdgesOnCellSrc->ptr1D(), verticesOnCellSrc->ptr2D(), cellsOnVertexSrc->ptr2D(),
                                      latCellSrc->ptr1D(), lonCellSrc->ptr1D(), latVertexSrc->ptr1D(), lonVertexSrc->ptr1D(), zmidSrc->ptr2D(),
                                      latCellDst->ptr1D(), lonCellDst->ptr1D(), zgridDst->ptr2D());
	stop_timer(0, &secs, &nsecs);
	printf("Time to create cellLevelMap : %i.%9.9i\n", secs, nsecs);


	//
	// Compute weights for remapping cell-based fields on layers
	//
	start_timer(0);
	cellLayerMap = new RemapperCell();
	cellLayerMap->computeWeightsCell(latCellDst->dimSize("nCells"), zmidSrc->dimSize("nVertLevels"), zmidDst->dimSize("nVertLevels"), 3,
                                      nEdgesOnCellSrc->ptr1D(), verticesOnCellSrc->ptr2D(), cellsOnVertexSrc->ptr2D(),
                                      latCellSrc->ptr1D(), lonCellSrc->ptr1D(), latVertexSrc->ptr1D(), lonVertexSrc->ptr1D(), zmidSrc->ptr2D(),
                                      latCellDst->ptr1D(), lonCellDst->ptr1D(), zmidDst->ptr2D());
	stop_timer(0, &secs, &nsecs);
	printf("Time to create cellLayerMap : %i.%9.9i\n", secs, nsecs);


	//
	// Compute weights for remapping cell-based fields to edges (used for interpolating rho to rho_edge)
	//
	start_timer(0);
	cellToEdgeMap = new RemapperCell();
	cellToEdgeMap->computeWeightsCell(latEdgeDst->dimSize("nEdges"), zedgeSrc->dimSize("nVertLevels"), zedgeDst->dimSize("nVertLevels"), 3,
                                      nEdgesOnCellSrc->ptr1D(), verticesOnCellSrc->ptr2D(), cellsOnVertexSrc->ptr2D(),
                                      latCellSrc->ptr1D(), lonCellSrc->ptr1D(), latVertexSrc->ptr1D(), lonVertexSrc->ptr1D(), zmidSrc->ptr2D(),
                                      latEdgeDst->ptr1D(), lonEdgeDst->ptr1D(), zedgeDst->ptr2D());
	stop_timer(0, &secs, &nsecs);
	printf("Time to create cellToEdgeMap : %i.%9.9i\n", secs, nsecs);


	//
	// Compute weights for remapping edge-based fields on layers (only used for u and v right now)
	//
	start_timer(0);
	edgeMap = new RemapperEdge();
	edgeMap->computeWeightsEdge(latCellSrc->dimSize("nCells"), latEdgeDst->dimSize("nEdges"),
                                      zedgeSrc->dimSize("nVertLevels"), zedgeDst->dimSize("nVertLevels"),
                                      nEdgesOnCellSrc->ptr1D(), cellsOnCellSrc->ptr2D(), edgesOnCellSrc->ptr2D(),
                                      latCellSrc->ptr1D(), lonCellSrc->ptr1D(), latEdgeSrc->ptr1D(), lonEdgeSrc->ptr1D(), zedgeSrc->ptr2D(),
                                      latCellDst->ptr1D(), lonCellDst->ptr1D(), latEdgeDst->ptr1D(), lonEdgeDst->ptr1D(), zedgeDst->ptr2D());
	stop_timer(0, &secs, &nsecs);
	printf("Time to create edgeMap : %i.%9.9i\n", secs, nsecs);


	//
	// Time-dependent processing for all global input times
	//
	for (itime=argv_idx; itime<argc; itime++) {
		globalFieldFile = argv[itime];


		//
		// Allocate and read global input fields
		//
		start_timer(0);
		try {
			xtime = new NCField<char>(globalFieldFile, "xtime");
		}
		catch (int e) {
			std::cerr << "Error reading xtime field from " << globalFieldFile << std::endl;
			return 1;
		}
		uSrc = new NCField<float>(globalFieldFile, "u");
		vSrc = new NCField<float>("v", 3, "Time", (size_t)1, "nEdges", uSrc->dimSize("nEdges"), "nVertLevels", uSrc->dimSize("nVertLevels"));
		theta_mSrc = new NCField<float>(globalFieldFile, "theta_m");
		rho_zzSrc = new NCField<float>(globalFieldFile, "rho_zz");
		wSrc = new NCField<float>(globalFieldFile, "w");
		stop_timer(0, &secs, &nsecs);
		printf("Time to read time-dependent fields from %s : %i.%9.9i\n", globalFieldFile, secs, nsecs);


		//
		// Allocate fields for interpolated regional fields
		//
		uDst = new NCField<float>("lbc_u", 3, "Time", (size_t)1, "nEdges", angleEdgeDst->dimSize("nEdges"), "nVertLevels", zmidDst->dimSize("nVertLevels"));
		vDst = new NCField<float>("lbc_v", 3, "Time", (size_t)1, "nEdges", angleEdgeDst->dimSize("nEdges"), "nVertLevels", zmidDst->dimSize("nVertLevels"));
		theta_mDst = new NCField<float>("lbc_theta_m", 3, "Time", (size_t)1, "nCells", zmidDst->dimSize("nCells"), "nVertLevels", zmidDst->dimSize("nVertLevels"));
		rho_zzDst = new NCField<float>("lbc_rho_zz", 3, "Time", (size_t)1, "nCells", zmidDst->dimSize("nCells"), "nVertLevels", zmidDst->dimSize("nVertLevels"));
		rho_edgeDst = new NCField<float>("lbc_rho_edge", 3, "Time", (size_t)1, "nEdges", angleEdgeDst->dimSize("nEdges"), "nVertLevels", zmidDst->dimSize("nVertLevels"));
		wDst = new NCField<float>("lbc_w", 3, "Time", (size_t)1, "nCells", zmidDst->dimSize("nCells"), "nVertLevelsP1", zgridDst->dimSize("nVertLevelsP1"));


		uSrcArr = uSrc->ptr3D();
		vSrcArr = vSrc->ptr3D();
		nEdgesOnEdgeSrcArr = nEdgesOnEdgeSrc->ptr1D();
		edgesOnEdgeSrcArr = edgesOnEdgeSrc->ptr2D();
		weightsOnEdgeSrcArr = weightsOnEdgeSrc->ptr2D();
		angleEdgeSrcArr = angleEdgeSrc->ptr1D();


		//
		// Reconstruct the global v field, and rotate the {u,v} vector field so that
		// u is the zonal wind component and v is the meridional wind component
		//
		start_timer(0);
		reconstruct_v(uSrc->dimSize("nEdges"), uSrc->dimSize("nVertLevels"), nEdgesOnEdgeSrcArr, edgesOnEdgeSrcArr, weightsOnEdgeSrcArr, uSrcArr[0], vSrcArr[0]);
		rotate_winds(uSrc->dimSize("nEdges"), uSrc->dimSize("nVertLevels"), angleEdgeSrcArr, uSrcArr[0], vSrcArr[0], 1);
		stop_timer(0, &secs, &nsecs);
		printf("Time to reconstruct v and rotate winds : %i.%9.9i\n", secs, nsecs);

		zzSrcArr = zzSrc->ptr2D();
		zzDstArr = zzDst->ptr2D();

		rho_zzSrcArr = rho_zzSrc->ptr3D();
		rho_zzDstArr = rho_zzDst->ptr3D();


		//
		// Uncouple global rho_zz field
		//
		start_timer(0);
		uncouple(zzSrc->dimSize("nVertLevels"), zzSrc->dimSize("nCells"), rho_zzSrcArr[0], zzSrcArr);
		stop_timer(0, &secs, &nsecs);
		printf("Time to uncouple rho_zz : %i.%9.9i\n", secs, nsecs);


		//
		// Set up name of regional output file as lbc.yyyy-mm-dd_hh.nc
		//
		xtimeArr = xtime->ptr2D();
		snprintf(date, (size_t)14, "%s", xtimeArr[0]);
		snprintf(regionalFieldFile, (size_t)64, "lbc.%s.nc", date);


		//
		// Create output file and define fields in it
		//
		stat = nc_create(regionalFieldFile, NC_64BIT_OFFSET, &ncid);

		stat = xtime->defineInFile(ncid);
		stat = uDst->defineInFile(ncid);
		stat = theta_mDst->defineInFile(ncid);
		stat = rho_zzDst->defineInFile(ncid);
		stat = rho_edgeDst->defineInFile(ncid);
		stat = wDst->defineInFile(ncid);

		//
		// Look for scalars to process (qv, qc, qr, etc.)
		//
		for (int i=0; i<NUM_SCALARS; i++) {
			try {
				qxSrc[i] = new NCField<float>(globalFieldFile, qxNames[i]);
				std::cout << "found " << qxNames[i] << " in " << globalFieldFile << std::endl;
				snprintf(qxLbcName, (size_t)64, "lbc_%s", qxNames[i]);
				qxDst[i] = new NCField<float>(qxLbcName, 3, "Time", (size_t)1, "nCells", zmidDst->dimSize("nCells"), "nVertLevels", zmidDst->dimSize("nVertLevels"));
				qxDst[i]->remapFrom(*qxSrc[i], *cellLayerMap);
				stat = qxDst[i]->defineInFile(ncid);
				delete qxSrc[i];
			}
			catch (int e) {
				std::cout << qxNames[i] << " not found in " << globalFieldFile << std::endl;
				qxDst[i] = new NCField<float>();
			}
		}

		stat = nc_enddef(ncid);

		//
		// Interpolate the zonal and meridional winds, and rotate the wind vector field so that u is the normal component
		//
		start_timer(0);
		uDstArr = uDst->ptr3D();
		vDstArr = vDst->ptr3D();
		angleEdgeDstArr = angleEdgeDst->ptr1D();
		uDst->remapFrom(*uSrc, *edgeMap);
		vDst->remapFrom(*vSrc, *edgeMap);
		rotate_winds(uDst->dimSize("nEdges"), uDst->dimSize("nVertLevels"), angleEdgeDstArr, uDstArr[0], vDstArr[0], 0);


		//
		// Interpolate scalar fields
		//
		theta_mDst->remapFrom(*theta_mSrc, *cellLayerMap);
		rho_zzDst->remapFrom(*rho_zzSrc, *cellLayerMap);
		rho_edgeDst->remapFrom(*rho_zzSrc, *cellToEdgeMap);
		wDst->remapFrom(*wSrc, *cellLevelMap);
		stop_timer(0, &secs, &nsecs);
		printf("Time to remap fields : %i.%9.9i\n", secs, nsecs);

		start_timer(0);


		//
		// Couple rho_zz
		//
		couple(zzDst->dimSize("nVertLevels"), zzDst->dimSize("nCells"), rho_zzDstArr[0], zzDstArr);


		//
		// Couple rho_edge
		//
		zzEdgeDstArr = zzEdgeDst->ptr2D();
		rho_edgeDstArr = rho_edgeDst->ptr3D();
		avg_cell_to_edge(zzEdgeDst->dimSize("nEdges"), zzEdgeDst->dimSize("nVertLevels"), cellsOnEdgeDst->ptr2D(), zzDstArr, zzEdgeDstArr);
		couple(zzEdgeDst->dimSize("nVertLevels"), zzEdgeDst->dimSize("nEdges"), rho_edgeDstArr[0], zzEdgeDstArr);
		stop_timer(0, &secs, &nsecs);
		printf("Time to couple rho_zz and rho_edge : %i.%9.9i\n", secs, nsecs);


		//
		// Write interpolated regional fields to output file
		//
		start_timer(0);
		stat = xtime->writeToFile(ncid);
		stat = uDst->writeToFile(ncid);
		stat = theta_mDst->writeToFile(ncid);
		stat = rho_zzDst->writeToFile(ncid);
		stat = rho_edgeDst->writeToFile(ncid);
		stat = wDst->writeToFile(ncid);
		stop_timer(0, &secs, &nsecs);

		for (int i=0; i<NUM_SCALARS; i++) {
			if (qxDst[i]->valid()) {
				stat = qxDst[i]->writeToFile(ncid);
			}
			delete(qxDst[i]);
		}

		printf("Time to write output fields : %i.%9.9i\n", secs, nsecs);

		stat = nc_close(ncid);


		delete xtime;
		delete uSrc;
		delete vSrc;
		delete theta_mSrc;
		delete rho_zzSrc;
		delete wSrc;

		delete uDst;
		delete vDst;
		delete theta_mDst;
		delete rho_zzDst;
		delete rho_edgeDst;
		delete wDst;
	}
	

	delete cellsOnEdgeDst;
	delete latCellDst;
	delete lonCellDst;
	delete latEdgeDst;
	delete lonEdgeDst;
	delete angleEdgeDst;
	delete zgridDst;
	delete zedgeDst;
	delete zmidDst;
	delete latCellSrc;
	delete lonCellSrc;
	delete latEdgeSrc;
	delete lonEdgeSrc;
	delete latVertexSrc;
	delete lonVertexSrc;
	delete nEdgesOnCellSrc;
	delete nEdgesOnEdgeSrc;
	delete weightsOnEdgeSrc;
	delete edgesOnEdgeSrc;
	delete angleEdgeSrc;
	delete cellsOnCellSrc;
	delete verticesOnCellSrc;
	delete cellsOnVertexSrc;
	delete cellsOnEdgeSrc;
	delete edgesOnCellSrc;
	delete zgridSrc;
	delete zedgeSrc;
	delete zmidSrc;
	delete zzSrc;
	delete zzDst;
	delete zzEdgeDst;
	delete cellLayerMap;
	delete cellLevelMap;
	delete edgeMap;
	delete cellToEdgeMap;

	return 0;
}
