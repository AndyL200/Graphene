#pragma once
#include <ft2build.h>
#include FT_FREETYPE_H

#include "Core.h"
#include "ShaderHead.h"
#include "pchgrfy.h"

//using method from learnopengl.com

struct Character
{
	unsigned int TextureID;
	glm::ivec2 Size;
	glm::ivec2 Bearing;
	unsigned int Advance;
};

Character* loadFont(const char* fontFileName, unsigned int fontSize)
{
	FT_Library ft;
	if (FT_Init_FreeType(&ft))
	{
		std::cout << "Error: could not init";
		return nullptr;
	}

	FT_Face face;
	if (FT_New_Face(ft, "fonts/" + *fontFileName, 0, &face))
	{
		std::cout << "ERROR:FREETYPE: Failed to load font";
		return nullptr;
	}

	//expand later for character sets from other languages
	FT_Set_Pixel_Sizes(face, 0, fontSize);
	Character characters[128];

	for (unsigned char c = 0; c < 127; c++)
	{
		if (FT_Load_Char(face, c, FT_LOAD_RENDER))
		{
			std::cout << "ERROR: Failed to load glyph" << std::endl;
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
	//REQUIRES ALL TEXTURES HAVE 4-BYTE ALIGNMENT
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	FT_Done_Face(face);
	FT_Done_FreeType(ft);

	return characters;
}

void renderText(Shader &s, std::string text, const char* font, float x, float y, float scale, glm::vec3 color) {

	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	
	VAO VAO1;
	VBO VBO1;

	s.Activate();
	glUniform3f(glGetUniformLocation(s.ID, "textColor"), color.x, color.y, color.z);
	glActiveTexture(GL_TEXTURE0);

	VAO1.Bind();


	Character* font_chars = loadFont(font, (int)scale);

	std::string::const_iterator c;
	for (c = text.begin(); c != text.end(); c++)
	{
		Character ch = font_chars[*c];

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

		//render glyph texture over quad
		VAO1.Bind();
		VBO1.Bind();
		glBindTexture(GL_TEXTURE_2D, ch.TextureID);

		//updating content of VBO memory
		glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(vertices), vertices);
		glBindBuffer(GL_ARRAY_BUFFER, 0);

		glDrawArrays(GL_TRIANGLES, 0, 6);

		x += (ch.Advance >> 6) * scale;
		VBO1.Unbind();
	}
	VAO1.Unbind();

	glBindVertexArray(0);
	glBindTexture(GL_TEXTURE_2D, 0);
}