#ifndef _NCFIELDBASE_H
#define _NCFIELDBASE_H
class NCFieldBase
{
public:
	NCFieldBase();
	~NCFieldBase();
	virtual int readFromFile(const char *filename, const char *fieldname) = 0;
	virtual int readFromFile(int ncid, const char *fieldname) = 0;
	virtual int defineInFile(int ncid) = 0;
	virtual int writeToFile(int ncid) = 0;
};
#endif
