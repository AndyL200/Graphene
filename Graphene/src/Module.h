#ifndef GRAPHENE_MODULE
#define GRAPENE_MODULE
#include "Graphene.h"
#include <glad/glad.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

namespace Graphene {

	
	

	struct verts {
		//Modules should contain an cylic list of vertices to match XGrafix behavior with OpenGL
		GLfloat* coordinates;
		GLfloat* rgb;
		GLfloat* normal;
		verts(GLfloat* coor, GLfloat* rgb, GLfloat* norm)
		{
			coordinates = coor;
			rgb = rgb;
			normal = norm;
		};
	};

	GLfloat* raw_to_rgb(int raw) {
		GLfloat* r = malloc(sizeof(GLfloat) * 3);
		//deconstruct raw color data
	}

	ArrayList<verts> plot_structure_to_vertices(structureData s, dim d)
	{
		ArrayList<verts> vertices;
		int color = (s.fillFlag) ? s.fillColor : s.lineColor;
		if (d == TWOD)
		{
			unsigned int i = s.n;
			while (i > 0)
			{
				GLfloat* coor = malloc(sizeof(GLfloat)*2);
				coor[0] = s.points[i];
				coor[1] = s.points[i - 1];
				GLfloat* norm = malloc(sizeof(GLfloat) * 2);
				::glm::vec3 normal = ::glm::normalize(::glm::vec3(coor[0], coor[1], 0));
				norm[0] = normal.x;
				norm[1] = normal.y;
				struct verts v(coor, raw_to_rgb(color), norm); //glm function here)
				i -= 2;
				vertices.add(v);
			}
		}
		else {

			unsigned int i = s.n;
			while (i > 0)
			{
				GLfloat* coor;
				coor[0] = s.points[i];
				coor[1] = s.points[i - 1];
				coor[2] = s.points[i - 2];
				struct verts v(coor, raw_to_rgb(color)); //glm function here)
				i-=3;
				vertices.add(v);
			}
		}
		return vertices;
	}


	class Module
	{
	public:
		Module(FILE* filename, ArrayList<ArrayList<verts>> v);
		
		ArrayList<ArrayList<verts>> structures;
	};

	//extract plots from Module
	class Module_Plot : public Module
	{
	public:
		Module_Plot();

	};

}
#endif