#ifndef GRAPHENE_READING
#define GRAPHENE_READING

#include <memory>;
#include <vector>;
#include <cstdint>
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

int GHWrite()
{

}

	

#endif