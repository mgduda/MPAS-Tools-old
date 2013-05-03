#Assumes ${NETCDF} is defined to the root of the netcdf library

CXX=icpc
MPICC=mpicc
CFLAGS= -O3
EXECUTABLE= PointsToMpas.x

NCINCDIR=${NETCDF}/include
NCLIBDIR=${NETCDF}/lib
PNCINCDIR=${PNETCDF}/include
PNCLIBDIR=${PNETCDF}/lib
MPIINCDIR=/opt/ibmhpc/pe1209/mpich2/gnu/include64

all:
#	${MPICC} ${CFLAGS} -o ${EXECUTABLE} points-mpas.cpp -I${NCINCDIR} -I${PNCINCDIR} -L${NCLIBDIR} -L${PNCLIBDIR} -lpnetcdf -lnetcdf_c++ -lnetcdf
	${CXX} ${CFLAGS} -c points-mpas.cpp -I${NCINCDIR} -I${PNCINCDIR} -I${MPIINCDIR}
	${MPICC} ${CFLAGS} -o ${EXECUTABLE} points-mpas.o -L${NCLIBDIR} -L${PNCLIBDIR} -lpnetcdf -lnetcdf_c++ -lnetcdf

clean:
	rm -f grid.nc
	rm -f graph.info
	rm -f ${EXECUTABLE}

