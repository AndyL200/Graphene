#ifndef GRAPHENE_MODULE
#define GRAPENE_MODULE

#include <pchgrfy.h>
#include <glad/glad.h>

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

	GLfloat* raw_to_rgb(uint32_t raw) {
		char color[4] = { 0 };
		//deconstruct raw color data
		//int type is 4 bytes long
		//just bitshift the raw down?
		//based on endian-ness though

		ENDIAN e = getGlobalEndianess();
		if (e == BIG)
		{
			for (int byte = 3; byte >= 0; --byte) {
				for (int bit = 0; bit < 8; bit++)
				{
					color[byte] |= ((raw & 1) << bit);
					raw >>= 1;
				}
			}
		}
		else
		{
			for (int byte = 0; byte < 4; ++byte) {
				for (int bit = 0; bit < 8; ++bit)
				{
					color[byte] |= ((raw & 1) << bit);
					raw >>= 1;
				}
			}
		}
		GLfloat colour[4] = { 0 };
		colour[0] = static_cast<float>(color[0]);
		colour[1] = static_cast<float>(color[1]);
		colour[2] = static_cast<float>(color[2]);
		colour[3] = static_cast<float>(color[3]);
		return colour;
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
				GLfloat* coor = (GLfloat*)malloc(sizeof(GLfloat)*2);
				coor[0] = s.points[i];
				coor[1] = s.points[i - 1];
				GLfloat* norm = (GLfloat*)malloc(sizeof(GLfloat) * 2);
				::glm::vec3 normal = ::glm::normalize(::glm::vec3(coor[0], coor[1], 0));
				//no loss because values are normalized
				norm[0] = normal.x;
				norm[1] = normal.y;
				verts v(coor, raw_to_rgb(color), norm); //glm function here)
				i -= 2;
				vertices.add(&v);
			}
		}
		else {

			unsigned int i = s.n;
			while (i > 0)
			{
				GLfloat* coor = (GLfloat*)malloc(sizeof(GLfloat) * 3);
				coor[0] = s.points[i];
				coor[1] = s.points[i - 1];
				coor[2] = s.points[i - 2];
				GLfloat* norm = (GLfloat*)malloc(sizeof(GLfloat) * 3);
				::glm::vec3 normal = ::glm::normalize(::glm::vec3(coor[0], coor[1], 0));
				norm[0] = normal.x;
				norm[1] = normal.y;
				norm[2] = normal.z;
				verts v(coor, raw_to_rgb(color), norm); //glm function here)
				i-=3;
				vertices.add(&v);
			}
		}
		return vertices;
	}


	class Module
	{
	public:
		Module(const char* name, ArrayList<ArrayList<verts>> v);
		
		void writeToFile();
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