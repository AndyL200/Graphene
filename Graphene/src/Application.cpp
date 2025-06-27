#include "Application.h"


namespace Graphene {
	Application::Application() {
		
	}

	void onUpdate() {
		glfwPollEvents();
	}

	bool init(int argc, char* argv[], &t) {
		if (!glfwInit()) {
			std::cout << "INIT ERROR";
			return false;
		}
		return true;
	}

}