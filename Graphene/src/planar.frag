#version 330 core
in vec3 fragColor;
out vec4 final_color;

void main()
{
final_color = vec4(fragColor,1.0);
}