#ifndef GRAPHENE_READING



#define GRAPHENE_READING

#include <streambuf>
#include <vector>
#include <iostream>

namespace Graphene {

	enum class inputType
	{
		LINE_PLOT,				// compiling xgmovie was throwing error for unknown LINE_PLOT; JK 2019-01-14; 
		SCATTER_PLOT,
		VECTOR_PLOT,
		SURFACE_PLOT,
		SCATTER3D,
		IRREGULAR3D
	};


	class VectorStreamBuf : public std::streambuf {
	public:
		VectorStreamBuf(std::vector<char>& v);
		VectorStreamBuf() : std::streambuf() {}

		std::streambuf* setbuf(char* s, std::streamsize n) override;
	};
	class VectorStream : public std::iostream {
	public:

		VectorStream(std::vector<char>& v) : vbuf(v), std::iostream(&vbuf) {}

		VectorStream() : std::iostream(nullptr) {}

		VectorStreamBuf getVBuf();
	private:
		VectorStreamBuf vbuf;
	};

	//from buffer to ptr
	template <typename T>
	T* GHRead(VectorStream& in)
	{
		T* ptr;
		size_t size = sizeof(T);
		if (size <= sizeof((void*)0)) {
			return nullptr;
		}

		if (in.rdbuf()->in_avail() < size)
			return nullptr;

		ptr = reinterpret_cast<T*>(malloc(size));
		if (!ptr)
			return nullptr;
		//&ptr[0]
		in.read((char*)&ptr[0], size);
		return ptr;
	}
	//from ptr to buffer

	//can there be more than one of these objects in the buffer though?
	int GHWrite(void*& ptr, VectorStream& in);


	//so far only works for primitive types



	//to read in bytes
	std::vector<char> file_to_buffer(FILE* fp, size_t to_read, bool locator_change);


	inputType determineFType(FILE** fp, char* filename, int index);

}
#endif