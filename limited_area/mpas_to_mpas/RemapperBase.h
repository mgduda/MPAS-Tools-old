#ifndef _REMAPPERBASE_H
#define _REMAPPERBASE_H

#include <typeinfo>

class RemapperBase {
public:
	RemapperBase();
	~RemapperBase();
	virtual void remap(const std::type_info& t, int ndims, void *dst, void *src) = 0;
};
#endif
