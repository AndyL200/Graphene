#ifndef GRAPHENE_VBO
#define GRAPHENE_VBO

#include <glad/glad.h>

class VBO {
public:
	GLuint ID;
	VBO();

	VBO(GLfloat* vertices, GLsizeiptr size);

	void Bind();
	void Unbind();
	void Delete();


};

#endif