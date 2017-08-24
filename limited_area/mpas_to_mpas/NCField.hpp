#ifndef _NCFIELD_H
#define _NCFIELD_H

#include <iostream>
#include <exception>
#include <cstdarg>
#include <typeinfo>
#include <string.h>
#include "netcdf.h"
#include "NCFieldBase.h"
#include "RemapperBase.h"
#include "array_utils.hpp"

static const int MAXDIMS = 3;

template <class fieldType>
class NCField: virtual public NCFieldBase
{
public:
	NCField()
	{
		ndims = 0;
		for (int i=0; i<MAXDIMS; i++) {
			dimlens[i] = 0;
		}
		data = NULL;
		data2d = NULL;
		data3d = NULL;
		isValid = false;
		unlimited_dim = -1;
	}


	NCField(const char *filename, const char *fieldname)
	{
		int ncid;
		int stat;
	
		ndims = 0;
		for (int i=0; i<MAXDIMS; i++) {
			dimlens[i] = 0;
		}
		data = NULL;
		data2d = NULL;
		data3d = NULL;
		isValid = false;

		stat = nc_open(filename, NC_SHARE, &ncid);
		if (stat != NC_NOERR) {
			throw stat;
		}
	
		init(ncid, fieldname);

		stat = nc_close(ncid);
		if (stat != NC_NOERR) {
			throw stat;
		}

		isValid = true;
	}


	NCField(int ncid, const char *fieldname)
	{
		ndims = 0;
		for (int i=0; i<MAXDIMS; i++) {
			dimlens[i] = 0;
		}
		data = NULL;
		data2d = NULL;
		data3d = NULL;
		isValid = false;

		init(ncid, fieldname);

		isValid = true;
	}


	NCField(const char *fieldname_in, int ndims_in, ...)
	{
		va_list args;
		va_start(args, ndims_in);

		strncpy(varname, fieldname_in, NC_MAX_NAME);
		ndims = ndims_in;

		totsize = 1;
		for (int i=0; i<ndims; i++) {
			strncpy(dimnames[i], va_arg(args, char*), NC_MAX_NAME);
			dimlens[i] = va_arg(args, size_t);
			totsize *= dimlens[i];
		}
		
		if (std::is_same<float, fieldType>::value) {
			xtype = NC_FLOAT;
		}
		else if (std::is_same<int, fieldType>::value) {
			xtype = NC_INT;
		}
		else if (std::is_same<char, fieldType>::value) {
			xtype = NC_CHAR;
		}
		data = new fieldType[totsize];
		data2d = NULL;
		data3d = NULL;
		isValid = true;
		unlimited_dim = -1;

		va_end(args);
	}


	~NCField()
	{
		ndims = 0;
		isValid = false;
		for (int i=0; i<MAXDIMS; i++) {
			dimlens[i] = 0;
		}
		if (data != NULL) {
			delete[] data;
		}
		if (data2d != NULL) {
			deallocate_2d<fieldType>(data2d);
		}
		if (data3d != NULL) {
			deallocate_3d<fieldType>(data3d);
		}
	}


	bool valid()
	{
		return isValid;
	}


	inline size_t size(int idim=(-1))
	{
		if (idim == -1) {
			return totsize;
		}
		else {
			return dimlens[idim];
		}
	}


	inline size_t dimSize(const char *dim)
	{
		for (int i=0; i<ndims; i++) {
			if (strncmp(dim, dimnames[i], NC_MAX_NAME) == 0) {
				return dimlens[i];
			}
		}
		return (size_t)0;
	}


	inline fieldType at(int i)
	{
		if (i >= 0 && i < totsize && data != NULL) {
			return data[i];
		}
		else {
			throw i;
		}
	}


	inline fieldType at(int j, int i)
	{
		size_t offset;

		offset = i + dimlens[ndims-1] * j;

		if (offset >= 0 && offset < totsize && data != NULL) {
			return data[offset];
		}
		else {
			throw offset;
		}
	}


	inline fieldType at(int k, int j, int i)
	{
		size_t offset;

		offset = i + dimlens[ndims-1] * j + dimlens[ndims-2] * dimlens[ndims-1] * k;

		if (offset >= 0 && offset < totsize && data != NULL) {
			return data[offset];
		}
		else {
			throw offset;
		}
	}


	int readFromFile(const char *filename, const char *fieldname)
	{
		int ncid;
		int stat;

		stat = nc_open(filename, NC_SHARE, &ncid);
		if (stat != NC_NOERR) {
			throw stat;
		}

		init(ncid, fieldname);

		stat = nc_close(ncid);
		if (stat != NC_NOERR) {
			throw stat;
		}

		return 0;
	}


	int readFromFile(int ncid, const char *fieldname)
	{
		init(ncid, fieldname);

		return 0;
	}


	int defineInFile(int ncid)
	{
		int stat;
		int dimids[MAXDIMS];

		for (int i=0; i<ndims; i++) {
			stat = nc_inq_dimid(ncid, dimnames[i], &dimids[i]);
			if (stat != NC_NOERR) {
				if (i != unlimited_dim) {
					stat = nc_def_dim(ncid, dimnames[i], dimlens[i], &dimids[i]);
				}
				else {
					stat = nc_def_dim(ncid, dimnames[i], NC_UNLIMITED, &dimids[i]);
				}
			}
		}

		stat = nc_inq_varid(ncid, varname, &varid);
		if (stat != NC_NOERR) {
			stat = nc_def_var(ncid, varname, xtype, ndims, dimids, &varid);
		}
	}


	int writeToFile(int ncid)
	{
		int stat;
		size_t startp[MAXDIMS];

		if (unlimited_dim == -1) {
			stat = nc_put_var(ncid, varid, (void *)data);
		}
		else {
			for (int i=0; i<MAXDIMS; i++) {
				startp[i] = 0;
			}
			stat = nc_put_vara(ncid, varid, startp, dimlens, (void *)data);
		}
	}


	NCField<fieldType>& operator=(NCField<fieldType>& r)
	{
		if (this != &r) {
			if (this->ndims != r.ndims) {
				throw -100;
			}
			for (int i=0; i<this->ndims; i++) {
				if (this->dimlens[i] != r.dimlens[i]) {
					throw -100 + i;
				}
			}
		}
		memcpy(this->data, r.data, sizeof(fieldType) * this->totsize);
		return *this;
	}


	NCField<fieldType>& operator+=(NCField<fieldType>& r)
	{
		if (this->ndims != r.ndims) {
			throw -200;
		}
		for (int i=0; i<this->ndims; i++) {
			if (this->dimlens[i] != r.dimlens[i]) {
				throw -200 + i;
			}
		}

		for (size_t i=0; i<this->totsize; i++) {
			this->data[i] += r.data[i];
		}

		return *this;
	}


	fieldType * ptr1D()
	{
		if (data) {
			return data;
		}
		else {
			return NULL;
		}
	}


	fieldType ** ptr2D()
	{
		if (data2d) {
			return data2d;
		}
		else if (data && ndims == 2) {
			data2d = allocate_2d<fieldType>(dimlens[0], dimlens[1], data);
			return data2d;
		}
		else {
			return NULL;
		}
	}


	fieldType *** ptr3D()
	{
		if (data3d) {
			return data3d;
		}
		else if (data && ndims == 3) {
			data3d = allocate_3d<fieldType>(dimlens[0], dimlens[1], dimlens[2], data);
			return data3d;
		}
		else {
			return NULL;
		}
	}


	int rank()
	{
		return ndims;
	}


	void remapFrom(NCField<fieldType>& src, RemapperBase& map)
	{
		if (ndims == 2 && src.rank() == 2) {
			void *src2d;
			void *dst2d;

			src2d = src.ptr2D();
			dst2d = ptr2D();
			map.remap(typeid(fieldType), 2, dst2d, src2d);
		}
		else if (ndims == 3 && src.rank() == 3) {
			void *src3d;
			void *dst3d;

			src3d = src.ptr3D();
			dst3d = ptr3D();
			map.remap(typeid(fieldType), 3, dst3d, src3d);
		}
		else {
			throw "Either src and dst ranks are different or they are of unsupported rank";
		}
	}


private:
	void init(int ncid, const char *fieldname)
	{
		int stat;
		int varid;
		int *dimids;
		int nunlim;
		int iunlim;
	
		stat = nc_inq_varid(ncid, fieldname, &varid);
		if (stat != NC_NOERR) {
			throw stat;
		}
	
		stat = nc_inq_varndims(ncid, varid, &ndims);
		if (stat != NC_NOERR) {
			throw stat;
		}
	
		stat = nc_inq_vartype(ncid, varid, &xtype);
		if (stat != NC_NOERR) {
			throw stat;
		}
		if (std::is_same<float, fieldType>::value && xtype != NC_FLOAT) {
			throw -998;
		}
		else if (std::is_same<int, fieldType>::value && xtype != NC_INT) {
			throw -997;
		}
		else if (std::is_same<char, fieldType>::value && xtype != NC_CHAR) {
			throw -996;
		}
	
		dimids = new int[ndims];
	
		stat = nc_inq_vardimid(ncid, varid, dimids);
		if (stat != NC_NOERR) {
			throw stat;
		}
	
		unlimited_dim = -1;
		stat = nc_inq_unlimdims(ncid, &nunlim, NULL);
		if (stat != NC_NOERR || nunlim > 1) {
			throw -996;
		}
		if (nunlim == 1) {
			stat = nc_inq_unlimdims(ncid, NULL, &iunlim);
		}

		for (int i=0; i<ndims; i++) {
			stat = nc_inq_dimlen(ncid, dimids[i], &dimlens[i]);
			if (stat != NC_NOERR) {
				throw stat;
			}

			stat = nc_inq_dimname(ncid, dimids[i], dimnames[i]);
			if (stat != NC_NOERR) {
				throw stat;
			}

			if (iunlim == dimids[i]) {
				unlimited_dim = i;
			}
		}
		for (int i=ndims; i<MAXDIMS; i++) {
			dimlens[i] = 1;
		}
	
		delete[] dimids;
	
		totsize = 1;
		for (int i=0; i<ndims; i++) {
			totsize *= dimlens[i];
		}
	
		data = new fieldType[totsize];
	
		stat = nc_get_var(ncid, varid, (void *)data);
		if (stat != NC_NOERR) {
			throw stat;
		}

		stat = nc_inq_varname(ncid, varid, varname);
		if (stat != NC_NOERR) {
			throw stat;
		}

	}


	size_t dimlens[MAXDIMS];
	char dimnames[MAXDIMS][NC_MAX_NAME+1];
	char varname[NC_MAX_NAME+1];
	int ndims;
	int unlimited_dim;
	nc_type xtype;
	int varid;
	size_t totsize;
	fieldType *data;
	fieldType **data2d;
	fieldType ***data3d;
	bool isValid;
};
#endif
