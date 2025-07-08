#include "Window.h"

namespace Graphene {
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
				plotInit(window, type);
				return;
			}

			plotInit(window);
			return;
		}

		freeInit(window);
		return;

	}

	//theory, inline may help keep the context defined in Application
	void inline TWODWindow::plotInit(GLFWwindow* window, int type = (WINDOW_PLOT | WINDOW2D)) {

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

	void TWODWindow::freeInit(GLFWwindow* window, int type = (WINDOW_FREE | WINDOW2D)) {

	}

	void TWODWindow::SetVertices(vector<float>& v) { this->vertices = v; }

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