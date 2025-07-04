#ifndef GRAPHENE_CORE
#define GRAPHENE_CORE


#ifdef _WIN32
	#ifdef BUILD_DLL
		#define Graphene_API _declspec(dllexport)
	#else 
		#define Graphene_API _declspec(dllimport)
	#endif

#define grfy_mmap CreateFileMapping
#else
#define grfy_mmap(...) mmap(__VA_ARGS__)
	#error Graphene supports Windows
#endif

typedef enum { LITTLE, BIG } ENDIAN;

static ENDIAN getGlobalEndianess()
{
	int numtest = 1;
	char* test = (char*)&numtest;
	return (*test == (char)1) ? LITTLE : BIG;
}

#define BIT(x) 1 << x


#define SCALAR double
#define SCALAR_CHAR "double"
#define XG_SCALAR_DOUBLE

#endif