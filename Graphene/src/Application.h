#pragma once
#include <GLFW/glfw3.h>
#include <glad/glad.h>
#include "Graphene.h"
#include <Windows/Window.h>
#include "Module.h"



namespace Graphene {

	enum GrapheneColorDefaults : GL_RGBA
	{
		GrapheneGrey = {0.1f, 0.5f, 0.9f, 1.0f},

	};


	class Application {
	public:
		Application(int argc = 0, char* argv[] = nullptr) {}

		virtual void onUpdate() const = 0;

		virtual bool init(int argc, char* argv[]) const = 0;
	
		void GHSetupWindow(const char* Win_Title, int width, int height, uint16_t WinType, Dimensions* label_ptr, int PlotType);

		void GHSetupWindow(const char* Win_Title, int width, int height, int WinType, Dimensions* label_ptr);

		void GHSetupWindow(const char* Win_Title, int width, int height, int WinType, Dimensions label_ptr);

		void GHSet2D(const char* PlotType, const char* X_Label, const char* Y_Label,
			const char* State,
			int ulx, int uly, SCALAR X_Scale, SCALAR Y_Scale,
			int X_Auto_Rescale, int Y_Auto_Rescale, SCALAR X_Min,
			SCALAR X_Max, SCALAR Y_Min, SCALAR Y_Max);
		
		void GHModCreate();
		void GHModLoadSetup(const char* Title, uint32_t WinType, Dimensions* label_ptr, uint16_t PlotType);

		void GHStart();

		private:
			bool is_running;
			GLFWwindow* main;
			WindowProps mainProps;

	};
}