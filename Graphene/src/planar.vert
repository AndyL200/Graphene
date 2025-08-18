#version 330 core
layout (location = 0) in vec3 pos;
layout (location = 1) in vec3 vertCol;
out vec3 fragColor;

uniform mat4 projection;

void main()
{
gl_Position = projection * vec4(pos, 1.0);
fragColor = vertCol;
}