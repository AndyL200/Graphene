#ifndef GRAPHENE_APPLICATION
#define GRAPHENE_APPLICATION


#include "Module.h"

#include <GHWindow/GHWindow.h>


namespace Graphene {

	namespace GrapheneColorDefaults
	{
		struct GrapheneColor {
			float r, b, g;
		};

		constexpr GrapheneColor Grey = { 0.1f, 0.5f, 0.9f };
		constexpr GrapheneColor ambient_background_1 = {0.1f, 0.7f, 1.0f};


	};


	Graphene_API class Application {
	public:
		Application(int argc = 0, char* argv[] = nullptr);

		virtual void onUpdate() const;

		virtual bool init(int argc, char* argv[]) const;
	
		void GHSetupWindow(Window* Win_Title, uint16_t WinType, Dimensions* label_ptr, int PlotType);

		void GHSetupWindow(Window* Win_Title, uint16_t WinType, Dimensions* label_ptr);

		void GHSet2D(const char* PlotType, const char* X_Label, const char* Y_Label,
			const char* State,
			int ulx, int uly, SCALAR X_Scale, SCALAR Y_Scale,
			int X_Auto_Rescale, int Y_Auto_Rescale, SCALAR X_Min,
			SCALAR X_Max, SCALAR Y_Min, SCALAR Y_Max);
		void SetVsync();
		Module GHModCreatePlot(const char* theInputFile);
		void GHModLoadSetup(const char* Title, uint32_t WinType, Dimensions* label_ptr, uint16_t PlotType);

		void GHStart();
		void GHClose();

		private:
			bool is_running;
			GLFWwindow* main;
			WindowProps mainProps;

	};
}


#endif