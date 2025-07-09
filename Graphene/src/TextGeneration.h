#include <ft2build.h>
#include FT_FREETYPE_H

#include <iostream>
#include <glm/glm.hpp>
#include <glad/glad.h>

struct Character
{
	unsigned int TextureID;
	glm::ivec2 Size;
	glm::ivec2 Bearing;
	unsigned int Advance;
};
FT_Face* loadFace(const char* fontFileName) {
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

	return &face;
}

Character* getCharacters(FT_Face face, unsigned int fontSize)
{
	FT_Set_Pixel_Sizes(face, 0, fontSize);
	Character characters[96];

	for (unsigned char c = 32; c < 127; c++)
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
			characters[c - 32] = tmpc;
	}

}