

#ifndef GRAPHENE_CORE
#define GRAPHENE_CORE


#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#ifdef _WIN32
	#ifdef BUILD_DLL
		#define GRAPHENE_API _declspec(dllexport)
	#else 
		#define GRAPHENE_API _declspec(dllimport)
	#endif

#define GH_PLATFORM_WINDOWS

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



enum dim {
	TWOD,
	THREED
};



#define SCALAR double
#define SCALAR_CHAR "double"
#define XG_SCALAR_DOUBLE

#endif