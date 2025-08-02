#ifndef GRAPHENE_MODULE
#define GRAPENE_MODULE

#include <pchgrfy.h>
#include "Plots.h"

namespace Graphene {


	

	struct verts {
		//Modules should contain an cylic list of vertices to match XGrafix behavior with OpenGL
		GLfloat* coordinates;  
		GLfloat* rgb;
		GLfloat* normal;
		verts(GLfloat* coor, GLfloat* rgb, GLfloat* norm)
		{
			this->coordinates = coor;
			this->rgb = rgb;
			this->normal = norm;
		};
	};
	//change to vector
	GLfloat* raw_to_rgb(uint32_t raw);

	ArrayList<verts> plot_structure_to_vertices(structureData s, dim d);


	class Module
	{
	public:
		Module(const char* name, ArrayList<ArrayList<verts>>& v);
		
		void writeToFile();
		ArrayList<ArrayList<verts>> structures;
	};

	//extract plots from Module
	/*class Module_Plot : public Module
	{
	public:
		Module_Plot();

	};*/

}
#endif