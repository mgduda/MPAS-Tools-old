void couple(int nVertLevels, int nCells, float **field, float **metric);
void uncouple(int nVertLevels, int nCells, float **field, float **metric);
void reconstruct_v(int nEdges, int nVertLevels, int *nEdgesOnEdge, int **edgesOnEdge, float **weightsOnEdge, float **u, float **v);
void rotate_winds(int nEdges, int nVertLevels, float *angleEdge, float **u, float **v, int toEarth);
void avg_to_midpoint(int nCells, int nLevels, float **levels, float **layers);
void avg_cell_to_edge(int nEdges, int nLevels, int **cellsOnEdge, float **cellfield, float **edgefield);
