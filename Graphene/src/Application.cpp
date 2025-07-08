#include "Application.h"


namespace Graphene {

	Application::Application(int argc = 0, char* argv[] = nullptr) {
		glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
		glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
		glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

		init(argc, argv);
		gladLoadGL();

		//read logic here



		
	}

	void Application::onUpdate() const {
		glfwPollEvents();
	}

	bool Application::init(int argc, char* argv[]) const {
		if (!glfwInit()) {
			std::cout << "INIT ERROR";
			return false;
		}
		return true;
	}
	bool Application::VSync() {
		glfwSwapInterval(1);
	}


	//assumes plot
	void Application::GHSetupWindow(Window* Win_Title, uint16_t WinType, Dimensions* label_ptr, int PlotType) {
		int width = Win_Title->GetWidth();
		int height = Win_Title->GetHeight();

		main = glfwCreateWindow(width, height, Win_Title->GetTitle(), NULL, NULL);

		glfwMakeContextCurrent(main);
		Win_Title->plotInit(main, WinType);


	}
	//assumes not plot
	void Application::GHSetupWindow(Window* Win_Title, int width, int height, uint16_t WinType, Dimensions* label_ptr) {
		int width = Win_Title->GetWidth();
		int height = Win_Title->GetHeight();

		main = glfwCreateWindow(width, height, Win_Title->GetTitle(), NULL, NULL);

		glfwMakeContextCurrent(main);
		Win_Title->freeInit();

	}

	void Application::GHSet2D(const char* PlotType, const char* X_Label, const char* Y_Label,
		const char* State,
		int ulx, int uly, SCALAR X_Scale, SCALAR Y_Scale,
		int X_Auto_Rescale, int Y_Auto_Rescale, SCALAR X_Min,
		SCALAR X_Max, SCALAR Y_Min, SCALAR Y_Max)
	{
		int             plottype;

		Dimensions* label_ptr;

		//label_ptr = SetupLabelStruct();
		label_ptr->Y_Label = strdup(Y_Label);
		label_ptr->Y_Min = Y_Min;
		label_ptr->Y_Max = Y_Max;
		label_ptr->Y_Scale = Y_Scale;
		label_ptr->Y_Auto_Rescale = Y_Auto_Rescale;

		label_ptr->X_Label = strdup(X_Label);
		label_ptr->X_Min = X_Min;
		label_ptr->X_Max = X_Max;
		label_ptr->X_Scale = X_Scale;
		label_ptr->X_Auto_Rescale = X_Auto_Rescale;

		if (!strcmp(PlotType, "linlin"))      plottype = LIN_LIN;
		else if (!strcmp(PlotType, "linlog")) plottype = LIN_LOG;
		else if (!strcmp(PlotType, "loglog")) plottype = LOG_LOG;
		else if (!strcmp(PlotType, "loglin")) plottype = LOG_LIN;
		else {
			printf("Unrecognized plot string '%s' for Window '%s'.\n", PlotType, Y_Label);
			exit(-1);
		}
		//TODO(ANDREW)PUT A FLAG HERE FOR LOAD VS CREATE
		GHSetupWindow(Y_Label, ulx, uly, WINDOW2D, label_ptr, plottype);
	}


	Module Application::GHModCreatePlot(const char* theInputFile) {
	
		FILE* inputfile;
		FILE* tmp;
		ArrayList<ArrayList<verts>> v;
		//convert plots to vertices?
		ArrayList<plot>* thePlots = new ArrayList<plot>;
		char line[512];
		char filename[512];
		char lastline[512];
		int xloc, yloc;
		//XInit(argc, argv, &theTime);
		if ((inputfile = fopen(theInputFile, "r")) == NULL) {
			fprintf(stderr, "Cannot open inputfile.  Exiting.");
			exit(1);
		}
		//theTime = 0;
		//maxTime = 0;
		//minTime = 1e9;
		fgets(line, 511, inputfile);
#ifdef XG_SCALAR_DOUBLE
		//sscanf(line, "%lf %d", &thetimestep, &MAX_LIN);
#else
		//sscanf(line, "%f %d", &thetimestep, &MAX_LIN);
#endif
		if (MAX_LIN <= 0) MAX_LIN = 12000;

		while (!feof(inputfile)) {
			fgets(line, 511, inputfile);
			sscanf(line, "%s %d %d", filename, &xloc, &yloc);
			//	 sscanf(line,"%d %d",&xloc,&yloc);
			// 

			xloc = std::max(xloc, 1); //xloc = min(700, xloc);
			yloc = std::max(yloc, 1); //yloc = min(700, yloc);
			if (!strcmp(filename, "END")) break;
			inputType t = openFile(&tmp, filename, 0);
			VectorStream inV = file_to_buffer(tmp, 0, true);
			switch (t)
			{
			
			case LINE_PLOT:
				LinePlot curPlot(inV, xloc, yloc);
				ArrayList<structureData> temp = curPlot.getStructures();
				for (int i = 0; i < temp.nItems; i++) {
					v.add(&plot_structure_to_vertices(temp[i], TWOD));
				}
				break;
			case SCATTER_PLOT:
				ScatterPlot curPlot(inV, xloc, yloc));
				ArrayList<structureData> temp = curPlot.getStructures();
				for (int i = 0; i < temp.nItems; i++) {
					v.add(&plot_structure_to_vertices(temp[i], TWOD));
				}
				break;
			case SURFACE_PLOT:
				SurfacePlot curPlot(inV, xloc, yloc));
				ArrayList<structureData> temp = curPlot.getStructures();
				for (int i = 0; i < temp.nItems; i++) {
					v.add(&plot_structure_to_vertices(temp[i], TWOD));
				}
				break;
			case VECTOR_PLOT:
				VectorPlot curPlot(inV, xloc, yloc));
				ArrayList<structureData> temp = curPlot.getStructures();
				for (int i = 0; i < temp.nItems; i++) {
					v.add(&plot_structure_to_vertices(temp[i], TWOD));
				}
				break;
			default:
				fprintf(stderr, "Unsupported plot type in file %s\n", filename);
			}
			fclose(tmp);
		}
		//create file using filename
		Module m(, v);
		return m;
	};
	void Application::GHModLoadSetup(const char* Title, uint32_t WinType, Dimensions* label_ptr, uint16_t PlotType)
	{
		char type;
		char tmpComm[501];

		glfwSetWindowTitle(main, Title);


		/*if (numberOfWindows > sizeOfWindowArray - 1) {
			ReallocateWindows();
		}*/

		//theWindowArray[numberOfWindows++] = theNewWindow;

		if (WinType & WINDOW3D) {
			THREEDWindow win(mainProps, PlotType);
			ColorCode_On(theNewWindow);
			Grid_On(theNewWindow);
			win->paint_function = Paint_ThreeD_Window;
			win->print_function = PostScript_ThreeD_Window;
			win->ascii_print_function = Ascii_ThreeD_Window;
			win->xgrafix_print_function = Bin_ThreeD_Window;
			if (PlotType & 0b100)
				Set_X_Log(theNewWindow);
			if (PlotType & 0b010)
				Set_Y_Log(theNewWindow);
			if (PlotType & 0b001)
				Set_Z_Log(theNewWindow);
		}
		else if (WinType == WINDOW_VECTOR_PLOT) {
			VectorWindow win(mainProps, PlotType);
			win->paint_function = Paint_Vector_Window;
			win->print_function = PostScript_Vector_Window;
			win->ascii_print_function = Ascii_Vector_Window;
			win->xgrafix_print_function = Bin_Vector_Window;
		}
		else {
			TWODWindow win(mainProps, PlotType);
			theNewWindow->paint_function = Paint_Window;
			theNewWindow->print_function = PostScriptOpenWindow;
			theNewWindow->ascii_print_function = Ascii_TwoD_Window;
			theNewWindow->xgrafix_print_function = Bin_TwoD_Window;
			if (PlotType & 0b10)
				Set_X_Log(theNewWindow);
			if (PlotType & 0b01)
				Set_Y_Log(theNewWindow);
		}

		if (label_ptr->Z_Auto_Rescale)
			Z_AutoRescale(theNewWindow);
		if (label_ptr->Y_Auto_Rescale)
			Y_AutoRescale(theNewWindow);
		if (label_ptr->X_Auto_Rescale)
			X_AutoRescale(theNewWindow);

		return;

}