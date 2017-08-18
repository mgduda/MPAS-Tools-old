#ifndef _ARRAY_UTILS_H
#define _ARRAY_UTILS_H
template <class T>
T ** allocate_2d(int dim1, int dim2, T *src)
{
	int i;
	T **arr;

	arr = (T **)malloc(sizeof(T *) * (size_t)dim1);
	for (i=0; i<dim1; i++) {
		arr[i] = src + (i * dim2);
	}

	return arr;
}


template <class T>
void deallocate_2d(T **arr)
{
	free(arr);
}


template <class T>
T *** allocate_3d(int dim1, int dim2, int dim3, T *src)
{
        size_t i, j;
        T ***arr3;
        T **arr2;

        arr3 = (T ***)malloc(sizeof(T **) * (size_t)dim1);
        arr2 = (T **)malloc(sizeof(T *) * (size_t)dim1 * (size_t)dim2);

        for (i=0; i<dim1; i++) {
                arr3[i] = arr2 + (i * dim2);
                for (j=0; j<dim2; j++) {
                        arr2[i*dim2+j] = src + (i*dim2*dim3) + (j * dim3);
                }
        }

        return arr3;
}


template <class T>
void deallocate_3d(T ***arr)
{
	free(arr[0]);
	free(arr);
}
#endif
