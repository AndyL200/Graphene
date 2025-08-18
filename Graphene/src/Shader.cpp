#include "Shader.h"

//#include <windows.h>

std::string get_file_contents(const char* filename)
{
	//wchar_t cwd[1000];
	//GetCurrentDirectory(1000, cwd);
	//std::wcout << cwd << std::endl;
	std::cout << filename << std::endl;

	std::ifstream in(filename, std::ios::binary);
	if (in)
	{

		std::string contents;
		in.seekg(0, std::ios::end);
		contents.resize(in.tellg());
		in.seekg(0, std::ios::beg);
		in.read(&contents[0], contents.size());
		in.close();
		return(contents);
	}
	std::cerr << "Error" << std::endl;
}

Shader::Shader(const char* vertexFile, const char* fragmentFile)
{
	std::string vertexCode = get_file_contents(vertexFile);

	std::string fragmentCode = get_file_contents(fragmentFile);

	const char* vertexSource = vertexCode.c_str();
	const char* fragmentSource = fragmentCode.c_str();

	GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);
	glShaderSource(vertexShader, 1, &vertexSource, NULL);
	glCompileShader(vertexShader);

	GLuint fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
	glShaderSource(fragmentShader, 1, &fragmentSource, NULL);
	glCompileShader(fragmentShader);

	ID = glCreateProgram();

	GLint success;
	GLchar infoLog[512];

	glAttachShader(ID, vertexShader);
	glAttachShader(ID, fragmentShader);


	glGetShaderiv(vertexShader, GL_COMPILE_STATUS, &success);
	if (!success) {
		glGetShaderInfoLog(vertexShader, 512, NULL, infoLog);
		std::cout << "Vertex Shader Error:\n" << infoLog << std::endl;
	}

	glGetShaderiv(fragmentShader, GL_COMPILE_STATUS, &success);
	if (!success) {
		glGetShaderInfoLog(fragmentShader, 512, NULL, infoLog);
		std::cout << "Fragment Shader Error:\n" << infoLog << std::endl;
	}


	glLinkProgram(ID);

	glGetProgramiv(ID, GL_LINK_STATUS, &success);
	if (!success) {
		glGetProgramInfoLog(ID, 512, NULL, infoLog);
		std::cout << "Shader Linking Error:\n" << infoLog << std::endl;
	}

	glDeleteShader(vertexShader);
	glDeleteShader(fragmentShader);
}

void Shader::Activate()
{
	glUseProgram(ID);
}
void Shader::Delete()
{
	glDeleteProgram(ID);
}

//void compileErrors()
//{
//
//}