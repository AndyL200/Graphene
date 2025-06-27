#pragma once

#ifdef _WIN32
	#ifdef BUILD_DLL
		#define Graphene_API _declspec(dllexport)
	#else 
		#define Graphene_API _declspec(dllimport)
	#endif
#else
	#error Graphene supports Windows
#endif

#define BIT(x) (1 << x);


#define SCALAR double
#define SCALAR_CHAR "double"
#define XG_SCALAR_DOUBLE