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
	}

