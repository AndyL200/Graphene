#ifndef GRAPHENE_READING



#define GRAPHENE_READING

#include <streambuf>
#include <vector>
#include <iostream>


class VectorStreamBuf : public std::streambuf {
public:
	VectorStreamBuf(std::vector<char>& v) { 
		char* begin = v.data();  
		size_t size = v.size();
		setbuf(begin, size);
	}
	VectorStreamBuf() : std::streambuf() {}

	std::streambuf* setbuf(char* s, std::streamsize n) override {
		setg(s, s, s + n);
		setp(s, s + n);
		return this;
	}
};
class VectorStream : public std::iostream {
public:

	VectorStream(std::vector<char>& v) : vbuf(v), std::iostream(&vbuf) {}

	VectorStream() : std::iostream(nullptr) {}

	VectorStreamBuf getVBuf()
	{
		return vbuf;
	}

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
	size_t size = sizeof(T);
	if (size <= sizeof((void*)0)) {
		return (T)0;
	}

	if (in.rdbuf()->in_avail() < size)
		return (T)0;

	ptr = reinterpret_cast<T*>(malloc(size));
	if (!ptr)
		return (T)0;

	in.read(ptr, size);
	return *ptr;
}

//to read in bytes
std::vector<char> file_to_buffer(FILE* fp, size_t to_read, bool locator_change)
{
	std::vector<char> v(to_read);
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
	return v;
}

inputType openFile(FILE** fp, char* filename, int index)
{
	char buf[100];
	//why is index even needed?
	if (index == 0 || index == 1) {
		if ((fopen_s(fp, filename, "r")) == NULL) {
			//sprintf(buf, "%s_%d", filename, index);
			//if ((*fp = fopen(buf, "r")) == NULL) {
				//fprintf(stderr, "Error, cannot find any file of the form %s\n", filename);
			exit(1);
		}
	}
	std::vector<char> sm = file_to_buffer(*fp, 1, true);
	VectorStream tmp(sm);
	inputType type = (inputType)GHRead<int>(tmp);
	return type;

}

	

#endif