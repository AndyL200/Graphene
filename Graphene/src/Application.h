#pragma once
#include <GLFW/glfw3.h>

namespace Graphene {
	class Application {
	public:
		Application() {}

		virtual void onUpdate() const = 0;

		virtual bool init() const = 0;
	private:
		bool is_running;
	};
}