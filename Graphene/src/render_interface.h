#pragma once

#include "ShaderHead.h"
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>


//at some point remove pragma onces

struct Glyph
{
	unsigned int textureID;
	unsigned int atlasWidth;
	unsigned int atlasHeight;
	glm::vec2 offset;
	glm::vec2 advance;
	glm::ivec2 txtCoords;
	glm::ivec2 size;
};

