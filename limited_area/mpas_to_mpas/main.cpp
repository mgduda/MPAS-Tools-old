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

	NCField<float> *latCellDst;
	NCField<float> *lonCellDst;
	NCField<float> *latEdgeDst;
	NCField<float> *lonEdgeDst;
	NCField<int> *cellsOnEdgeDst;
	NCField<char> *xtime;

	const char *globalMeshFile;
	const char *regionalMeshFile;
	const char *globalFieldFile;
	char regionalFieldFile[64];
	char **xtimeArr;

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
	NCField<float> *terSrc;
	NCField<float> *terDst;
	int *nEdgesOnEdgeSrcArr;
	RemapperCell *cellLevelMap;
	int secs, nsecs;


	if (argc < 3) {
		std::cerr << "\nUsage: " << argv[0] << " <global_static_file> <regional_static_file>\n";
		std::cerr << "       Terrain from the \"global\" file will be interpolated to the \"regional\" mesh\n\n";
		return 1;
	}

	//
	// Try to catch unintentional specification of options
	//
	if (argv[1][0] == '-') {
		std::cerr << "Unrecognized option: " << argv[1] << std::endl;
		std::cerr << "\nSupported options are: --use-reconstruct-winds\n";
		return 1;
	}

	globalMeshFile = argv[1];
	regionalMeshFile = argv[2];


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
	cellsOnEdgeDst = new NCField<int>(regionalMeshFile, "cellsOnEdge");
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
	cellsOnCellSrc = new NCField<int>(globalMeshFile, "cellsOnCell");
	verticesOnCellSrc = new NCField<int>(globalMeshFile, "verticesOnCell");
	cellsOnVertexSrc = new NCField<int>(globalMeshFile, "cellsOnVertex");
	cellsOnEdgeSrc = new NCField<int>(globalMeshFile, "cellsOnEdge");
	stop_timer(0, &secs, &nsecs);
	printf("Time to read mesh fields from %s : %i.%9.9i\n", globalMeshFile, secs, nsecs);


	//
	// Compute weights for remapping cell-based fields on levels (only used for w right now)
	//
	start_timer(0);
	cellLevelMap = new RemapperCell();
	cellLevelMap->computeWeightsCell(latCellDst->dimSize("nCells"), 0, 0, 3,
                                      nEdgesOnCellSrc->ptr1D(), verticesOnCellSrc->ptr2D(), cellsOnVertexSrc->ptr2D(),
                                      latCellSrc->ptr1D(), lonCellSrc->ptr1D(), latVertexSrc->ptr1D(), lonVertexSrc->ptr1D(), NULL,
                                      latCellDst->ptr1D(), lonCellDst->ptr1D(), NULL);
	stop_timer(0, &secs, &nsecs);
	printf("Time to create cellLevelMap : %i.%9.9i\n", secs, nsecs);


	//
	// Handle terrain field
	//
	{
		globalFieldFile = argv[1];


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

		try {
			terSrc = new NCField<float>(globalFieldFile, "ter");
		}
		catch (int e) {
			std::cerr << "Error reading ter field from " << globalFieldFile << std::endl;
			return 1;
		}

		stop_timer(0, &secs, &nsecs);
		printf("Time to read time-dependent fields from %s : %i.%9.9i\n", globalFieldFile, secs, nsecs);


		//
		// Allocate fields for interpolated regional fields
		//
		terDst = new NCField<float>("ter", 1, "nCells", latCellDst->dimSize("nCells"));


		//
		// Set up name of regional output file as terrain.nc
		//
		snprintf(regionalFieldFile, (size_t)64, "terrain.nc");


		//
		// Create output file and define fields in it
		//
		stat = nc_create(regionalFieldFile, NC_64BIT_OFFSET, &ncid);

		stat = terDst->defineInFile(ncid);

		stat = nc_enddef(ncid);

		//
		// Interpolate the zonal and meridional winds, and rotate the wind vector field so that u is the normal component
		//
		start_timer(0);
		terDst->remapFrom(*terSrc, *cellLevelMap);
		stop_timer(0, &secs, &nsecs);
		printf("Time to remap fields : %i.%9.9i\n", secs, nsecs);


		//
		// Write interpolated regional fields to output file
		//
		start_timer(0);
		stat = terDst->writeToFile(ncid);
		stop_timer(0, &secs, &nsecs);

		stat = nc_close(ncid);

		delete terDst;
	}
	

	delete terSrc;

	return 0;
}
