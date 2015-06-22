#Assumes ${NETCDF} is defined to the root of the netcdf library

CXX=mpicxx
MPICC=mpicxx
CFLAGS= -g
EXECUTABLE= PointsToMpas.x

NCINCDIR=${NETCDF}/include
NCLIBDIR=${NETCDF}/lib
PNCINCDIR=${PNETCDF}/include
PNCLIBDIR=${PNETCDF}/lib

all:
#	${MPICC} ${CFLAGS} -o ${EXECUTABLE} points-mpas.cpp -I${NCINCDIR} -I${PNCINCDIR} -L${NCLIBDIR} -L${PNCLIBDIR} -lpnetcdf -lnetcdf_c++ -lnetcdf
	${CXX} ${CFLAGS} -c points-mpas.cpp -I${NCINCDIR} -I${PNCINCDIR}
	${MPICC} ${CFLAGS} -o ${EXECUTABLE} points-mpas.o -L${NCLIBDIR} -L${PNCLIBDIR} -lpnetcdf -lnetcdf_c++ -lnetcdf

clean:
	rm -f grid.nc
	rm -f graph.info
	rm -f ${EXECUTABLE}

