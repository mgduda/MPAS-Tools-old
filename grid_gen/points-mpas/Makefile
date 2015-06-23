CXX=mpicxx
MPICC=mpicxx
CFLAGS= -g
EXECUTABLE= PointsToMpas.x

PNCINCDIR=${PNETCDF}/include
PNCLIBDIR=${PNETCDF}/lib

all:
	${CXX} ${CFLAGS} -c points-mpas.cpp -I${PNCINCDIR}
	${MPICC} ${CFLAGS} -o ${EXECUTABLE} points-mpas.o -L${PNCLIBDIR} -lpnetcdf

clean:
	rm -f ${EXECUTABLE} points-mpas.o

