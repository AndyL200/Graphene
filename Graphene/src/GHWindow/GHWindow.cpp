#include "GHWindow.h"

namespace Graphene {
	THREEDWindow::THREEDWindow(WindowProps props, int type) : Window(props) {
		this->type = type;
		//glfwWindowHint();
	}
	unsigned int THREEDWindow::GetWidth() {
		return props.Width;
	}
	unsigned int THREEDWindow::GetHeight() {
		return props.Height;
	}
	const char* THREEDWindow::GetTitle() {
		return props.Title.c_str();
	}

	/*GLFWwindow* THREEDWindow::Create()
	{

	}*/

	void THREEDWindow::plotInit(GLFWwindow* window) {
		VAO VAO1;

		gladLoadGL();
		int width, height;
		glfwGetWindowSize(window, &width, &height);
		glViewport(0, 0, width, height);


	}

	void THREEDWindow::freeInit(GLFWwindow* window) {

	}

	void THREEDWindow::OnUpdate() {
		glfwPollEvents();
	}
	TWODWindow::TWODWindow(WindowProps props, int type) : Window(props) {
		if (!(type & WINDOW2D)) {
			throw _EXCEPTION_("YOU CANT PARK THERE");
		}
		this->type = type;
		//glfwWindowHint();
	}

	GLFWwindow* TWODWindow::Create() {
		GLFWwindow* window = glfwCreateWindow(props.Width, props.Height, props.Title.c_str(), NULL, NULL);
		glfwMakeContextCurrent(window);
		if (type & WINDOW_PLOT)
		{
			if (type & (WINDOW_SCATTER_PLOT | WINDOW_VECTOR_PLOT | WINDOW_ONE_D_PLOT))
			{
				plotInit(window);
				return window;
			}

			plotInit(window);
			return window;
		}

		freeInit(window);
		return window;

	}

	//theory, inline may help keep the context defined in Application
	void inline TWODWindow::plotInit(GLFWwindow* window) {

		GLfloat verts[] = {
			// COORDINATES     // COLOR
				0.1f, 0.1f,		0.0f, 0.0f, 0.0f,
				1.0f, 0.1f,     1.0f, 0.0f, 0.0f,
				0.1f, 1.0f,     0.0f, 1.0f, 0.1f,
		};

		GLuint indices[] = {
			0, 1,
			0, 2
		};

		vertices = verts;

		VAO VAO1;


		VBO VBO1(vertices, sizeof(vertices));
		EBO EBO1(indices, sizeof(indices));
		VAO1.LinkAttrib(VBO1, 0, 2, GL_FLOAT, 5 * sizeof(float), (void*)0);
		VAO1.LinkAttrib(VBO1, 1, 3, GL_FLOAT, 5 * sizeof(float), (void*)(2 * sizeof(float)));

		//unbinding to prevent modification
		VBO1.Unbind();
		EBO1.Unbind();

		//typically add in your vertex shaders here


		glClear(GL_COLOR_BUFFER_BIT);

		glDrawElements(GL_LINES, sizeof(indices) / sizeof(int), GL_UNSIGNED_INT, indices);


		int width, height;
		glfwGetWindowSize(window, &width, &height);
		glViewport(0, 0, width, height);


	}

	void TWODWindow::freeInit(GLFWwindow* window) {

	}

	void TWODWindow::SetVertices(std::vector<float>& v) { this->vertices = v.data(); }

	GLfloat* TWODWindow::GetVertices() {
		return vertices;
	}

	unsigned int TWODWindow::GetWidth() {
		return props.Width;
	}
	unsigned int TWODWindow::GetHeight() {
		return props.Height;
	}
	const char* TWODWindow::GetTitle() {
		return props.Title.c_str();
	}

	void TWODWindow::OnUpdate() {
		glfwPollEvents();
	}
}

