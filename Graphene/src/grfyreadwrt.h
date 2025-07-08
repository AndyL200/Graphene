#ifndef GRAPHENE_READING
#define GRAPHENE_READING

#include <memory>;
#include <vector>;
#include <cstdint>;
#include "Core.h";
using std::vector;


class VectorStreamBuf : public std::streambuf {
public:
	VectorStreamBuf(vector<char>& v) { 
		char* begin = v.data();  
		size_t size = v.size();
		setg(begin, begin, begin + size);
		setp(begin, begin + size);
	}
};
class VectorStream : public std::iostream {
public:
	VectorStream(vector<char>& v) : std::iostream(&vbuf), vbuf(v) {};
private:
	VectorStreamBuf vbuf;
};

//from buffer to ptr
template <typename T>
T GHRead(VectorStream& in);
//from ptr to buffer

//can there be more than one of these objects in the buffer though?
int GHWrite(void*& ptr, VectorStream& in);

template <typename T>
T GHRead(VectorStream& in)
{
	T obj;
	T* ptr = &obj;
	VectorStreamBuf* buffer = in.rdbuf();
	size_t size = sizeof(T);
	if (size <= 0) {
		return -1;
	}

	if (buffer->in_avail() < size)
		return 1;

	ptr = malloc(size);
	if (!ptr)
		return 1;

	memcpy(ptr, buffer->gptr(), size);
	buffer->gbump(size);
	return obj;
}


int GHWrite(VectorStream& out, FILE* in)
{
	
}
//to read in bytes
VectorStream file_to_buffer(FILE* fp, size_t to_read, bool locator_change)
{
	vector<char> v(to_read);
	if (to_read == 0)
	{
		fseek(fp, 0, SEEK_END);
		long size = ftell(fp);
		rewind(fp);
		fread(v.data(), sizeof(char), size, fp);
	}
	else {
		fread(v.data(), sizeof(char), to_read, fp);
	}

	if (!locator_change)
	{
		rewind(fp);
	}
	return VectorStream(v);
}

inputType openFile(FILE** fp, char* filename, int index) {
	char buf[100];
	//why is index even needed?
	if (index == 0 || index == 1) {
		sprintf(buf, "%s", filename);
		if ((*fp = fopen(buf, "r")) == NULL) {
			//sprintf(buf, "%s_%d", filename, index);
			//if ((*fp = fopen(buf, "r")) == NULL) {
				//fprintf(stderr, "Error, cannot find any file of the form %s\n", filename);
			exit(1);
			}
		}
		VectorStream sm = file_to_buffer(*fp, 1, true);
		inputType type = (inputType)GHRead<int>(sm);
		return type;
	}
	

#endif