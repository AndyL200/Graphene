#ifndef GRAPHENE_APPLICATION
#define GRAPHENE_APPLICATION


#include "Module.h"
#include "GHWindow/GHWindow.h"




namespace Graphene {

	
		class GrapheneColor {
			//different formatting
		public:
			float r;
			float g;
			float b;
			GrapheneColor(float red, float green, float blue);
		};

		


	class GRAPHENE_API Application {
	public:
		Application(int argc = 0, char* argv[] = nullptr);
		//Application();

		void onUpdate();

		bool init(int argc, char* argv[]);
	
		void GHSetupWindow(Window* Win_Title, uint16_t WinType, Dimensions* label_ptr, int PlotType);

		void GHSetupWindow(Window* Win_Title, uint16_t WinType, Dimensions* label_ptr);

		void GHSet2D(const char* PlotType, const char* X_Label, const char* Y_Label,
			const char* State,
			int ulx, int uly, SCALAR X_Scale, SCALAR Y_Scale,
			int X_Auto_Rescale, int Y_Auto_Rescale, SCALAR X_Min,
			SCALAR X_Max, SCALAR Y_Min, SCALAR Y_Max);
		void SetVsync(bool b);
		Module GHModCreatePlot(const char* theInputFile);
		void GHModLoadSetup(const char* Title, uint32_t WinType, Dimensions* label_ptr, uint16_t PlotType);

		GLFWwindow* getWindow();

		void GHStart();
		void GHClose();

		private:
			bool is_running;
			GLFWwindow* main;
			WindowProps mainProps;

	};
}


#endif