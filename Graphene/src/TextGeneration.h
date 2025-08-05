#ifndef GRAPH_TEXT
#define GRAPH_TEXT

#include <ft2build.h>
#include FT_FREETYPE_H

#include "Core.h"
#include "Shader.h"
#include "VAO.h"
#include "VBO.h"
#include "EBO.h"
#include "pchgrfy.h"

	//using method from learnopengl.com

	struct Character
	{
		unsigned int TextureID;
		glm::ivec2 Size;
		glm::ivec2 Bearing;
		unsigned int Advance;
	};

	Character* loadFont(const char* fontFileName, unsigned int fontSize);

	void renderText(Shader& s, std::string text, const char* font, float x, float y, float scale, glm::vec3 color);

#endif