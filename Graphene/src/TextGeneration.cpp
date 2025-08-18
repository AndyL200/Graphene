#include "TextGeneration.h"


Character* loadFont(const char* fontFileName, unsigned int fontSize)
{
	FT_Library ft;
	if (FT_Init_FreeType(&ft))
	{
		std::cout << "Error: could not init";
		return nullptr;
	}

	FT_Face face;
	if (FT_New_Face(ft, fontFileName, 0, &face))
	{
		std::cout << "ERROR:FREETYPE: Failed to load font";
		return nullptr;
	}
	//REQUIRES ALL TEXTURES HAVE 4-BYTE ALIGNMENT
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	//expand later for character sets from other languages
	FT_Set_Pixel_Sizes(face, 0, fontSize);
	Character characters[128];

	for (unsigned char c = 0; c < 127; c++)
	{
		if (FT_Load_Char(face, c, FT_LOAD_RENDER))
		{
			std::cout << "ERROR: Failed to load glyph: " << c << std::endl;
			continue;
		}

		unsigned int texture;
		glGenTextures(1, &texture);
		glBindTexture(GL_TEXTURE_2D, texture);
		glTexImage2D(
			GL_TEXTURE_2D,
			0,
			GL_RED,
			face->glyph->bitmap.width,
			face->glyph->bitmap.rows,
			0,
			GL_RED,
			GL_UNSIGNED_BYTE,
			face->glyph->bitmap.buffer
		);

		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);


		Character tmpc = {
		texture,
		glm::ivec2(face->glyph->bitmap.width, face->glyph->bitmap.rows),
		glm::ivec2(face->glyph->bitmap_left, face->glyph->bitmap_top),
		face->glyph->advance.x
		};
		characters[c] = tmpc;
	}
	
	FT_Done_Face(face);
	FT_Done_FreeType(ft);

	return characters;
}


void renderText(Shader& s, std::string text, VAO v, Character* font_chars, float x, float y, float scale, glm::vec3 color)
{
	float copyX = x;
	s.Activate();
	
	glDisable(GL_DEPTH_TEST);
	
	glUniform3f(glGetUniformLocation(s.ID, "textColor"), color.x, color.y, color.z);
	glActiveTexture(GL_TEXTURE0);
	GLuint textLoc = glGetUniformLocation(s.ID, "text");

	VBO VBOTXT;
	v.LinkAttrib(VBOTXT, 0, 4, GL_FLOAT, 4 * sizeof(float), 0);
	glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 6 * 4, NULL, GL_DYNAMIC_DRAW);
	

	std::string::const_iterator c;
	for (c = text.begin(); c != text.end(); c++)
	{
		char care = *c;
		Character ch = font_chars[(int)*c];
		if (*c == '\n')
		{
			y -= ((ch.Size.y)) * 1.3 * scale;
			x = copyX;
		}
		else {

			//starts at the origin???
			float xpos = x + ch.Bearing.x * scale;
			float ypos = y - (ch.Size.y - ch.Bearing.y) * scale;

			float w = ch.Size.x * scale;
			float h = ch.Size.y * scale;

			float vertices[6][4] = {
				//Coords        //TexCoords
				{xpos, ypos + h,		0.0f, 0.0f},
				{xpos, ypos,			0.0f, 1.0f},
				{xpos + w, ypos,		1.0f, 1.0f},
				{xpos, ypos + h,		0.0f, 0.0f},
				{xpos + w, ypos,		1.0f, 1.0f},
				{xpos + w, ypos + h,	1.0f, 0.0f}
			};

			// render glyph texture over quad
			glBindTexture(GL_TEXTURE_2D, ch.TextureID);
			// update content of VBO memory
			VBOTXT.Bind();
			glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(vertices), vertices);
			// render quad
			glDrawArrays(GL_TRIANGLES, 0, 6);
			VBOTXT.Unbind();

			// now advance cursors for next glyph (note that advance is number of 1/64 pixels)
			x += (ch.Advance >> 6) * scale; // bitshift by 6 to get value in pixels (2^6 = 64)
		}
	}

	glBindVertexArray(0);
	glBindTexture(GL_TEXTURE_2D, 0);
}