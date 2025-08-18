#ifndef GH_SHADER
#define GH_SHADER

#include "pchgrfy.h"
#include <glad/glad.h>


class Shader
{
public:
	GLuint ID;
	Shader(const char* vertexFile, const char* fragmentFile);
	void Activate();
	void Delete();
	//void compileErrors(const char* type);
};

#endif