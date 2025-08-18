#include "Application.h"

#include <ft2build.h>
#include FT_FREETYPE_H

#include "render_interface.h"

std::vector<Glyph> loadFont(const char* fontPath, unsigned int fontSize);
void renderText(Shader& s, std::string text, VAO& v, std::vector<Glyph>& font_chars, float x, float y, float scale, glm::vec3 color);


std::vector<Glyph> loadFont(const char* fontPath, unsigned int fontSize)
{
	FT_Library ft;
	if (FT_Init_FreeType(&ft))
	{
		std::cout << "Error: could not init";
		return std::vector<Glyph>();
	}

	FT_Face face;
	if (FT_New_Face(ft, fontPath, 0, &face))
	{
		std::cout << "ERROR:FREETYPE: Failed to load font";
		return std::vector<Glyph>();
	}

	unsigned int texture;
	glGenTextures(1, (GLuint*)&texture);
	//REQUIRES ALL TEXTURES HAVE 4-BYTE ALIGNMENT
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	//expand later for character sets from other languages
	FT_Set_Pixel_Sizes(face, 0, fontSize);
	std::vector<Glyph> glyphs;

	int padding = 2;
	int row = 0;
	int col = padding;

	const int txtWidth = 512;
	char* textureBuffer = new char[txtWidth * txtWidth];

	for (FT_ULong glyphIdx = 32; glyphIdx < 127; ++glyphIdx)
	{
		FT_UInt glyphIndex = FT_Get_Char_Index(face, glyphIdx);
		FT_Load_Glyph(face, glyphIdx, FT_LOAD_DEFAULT);
		FT_Error error = FT_Render_Glyph(face->glyph, FT_RENDER_MODE_NORMAL);

		if (col + face->glyph->bitmap.width + padding >= 512)
		{
			col = padding;
			row += fontSize;
		}

		for (unsigned int y = 0; y < face->glyph->bitmap.rows; ++y)
		{
			for (unsigned int x = 0; x < face->glyph->bitmap.width; ++x)
			{
				textureBuffer[(row + y) * txtWidth + col + x] =
					face->glyph->bitmap.buffer[y * face->glyph->bitmap.pitch + x];
			}
		}
		Glyph g;
		g.textureID = texture;
		g.atlasWidth = 512;
		g.atlasHeight = 512;
		g.offset = glm::vec2((float)face->glyph->bitmap_left, (float)face->glyph->bitmap_top);
		g.advance = glm::vec2((float)(face->glyph->advance.x >> 6), (float)(face->glyph->advance.y >> 6));
		g.txtCoords = glm::vec2(col, row);
		g.size = glm::vec2(face->glyph->bitmap.width, face->glyph->bitmap.rows);

		col += face->glyph->bitmap.width + padding;
		glyphs.push_back(g);

		//std::cout << g.textureID << " " << g.atlasWidth << " " << g.offset.x << " " << g.offset.y << " " << g.advance.x << " " << g.advance.y << " " << g.txtCoords.x << " " << g.txtCoords.y << " " << g.size.x << " " << g.size.y << " " << std::endl;

	}



	//Upload OpenGL texture

	glActiveTexture(GL_TEXTURE1);
	glBindTexture(GL_TEXTURE_2D, texture);

	glTexImage2D(GL_TEXTURE_2D, 0, GL_R8, txtWidth, txtWidth, 0, GL_RED, GL_UNSIGNED_BYTE, (char*)textureBuffer);

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

	FT_Done_Face(face);
	FT_Done_FreeType(ft);


	return glyphs;
}


void renderText(Shader& s, std::string text, VAO& v, VBO& vb, std::vector<Glyph>& font_chars, float x, float y, float scale, glm::vec3 color)
{
	float copyX = x;
	//hardcoding texture 0 to be for text

	s.Activate();
	glUniform3f(glGetUniformLocation(s.ID, "textColor"), color.x, color.y, color.z);
	glUniform1i(glGetUniformLocation(s.ID, "text"), 1);
	glActiveTexture(GL_TEXTURE1);

	v.Bind();

	Glyph temp = font_chars[0];
	// render glyph texture over quad
	glBindTexture(GL_TEXTURE_2D, temp.textureID);


	


	//error zone
	VBO VBOTXT;
	glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 6 * 4, NULL, GL_DYNAMIC_DRAW);
	//first two aren't normalized next two are
	v.LinkAttrib(VBOTXT, 0, 2, GL_FLOAT, 4 * sizeof(float), 0);
	//error zone
	
	std::string::const_iterator c;
	for (c = text.begin(); c != text.end(); c++)
	{
		char care = *c;
		if (static_cast<int>(care) >= 32 && static_cast<int>(care) < 127)
		{
			Glyph g = font_chars[static_cast<int>(care) - 32];

			if (care == '\n')
			{
				y -= ((g.size.y)) * scale;
				x = copyX;
				continue;
			}


			//starts at the origin???
			float xpos = x + g.offset.x * scale;
			float ypos = y - (g.size.y - g.offset.y) * scale;

			float w = g.size.x * scale;
			float h = g.size.y * scale;

			float tx = g.txtCoords.x / g.atlasWidth;
			float ty = g.txtCoords.y / g.atlasHeight;
			float tw = g.size.x / g.atlasWidth;
			float th = g.size.y / g.atlasHeight;

			float vertices[6][4] = {
				//Coords        //TexCoords
				{xpos, ypos + h,		tx, ty},
				{xpos, ypos,			tx, ty + th},
				{xpos + w, ypos,		tx + tw, ty + th},
				{xpos, ypos + h,		tx, ty},
				{xpos + w, ypos,		tx + tw, ty + th},
				{xpos + w, ypos + h,	tx + tw, ty}
			};

			
			
			// update content of VBO memory
			//move out
			VBOTXT.Bind();
			glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(vertices), vertices);
			VBOTXT.Unbind();

			// render quad
			//std::cout << "issuing draw call" << std::endl;

			//error zone
			glDrawArrays(GL_TRIANGLES, 0, 6);
			//error zone
			GLenum err;

			while ((err = glGetError()) != GL_NO_ERROR)
			{
				std::cerr << "GL Error: " << err << std::endl;
			}
			//std::cout << "finished drawing" << std::endl;

			// now advance cursors for next glyph (note that advance is number of 1/64 pixels)
			x += g.advance.x * scale; // bitshift by 6 to get value in pixels (2^6 = 64)
		}

		glBindVertexArray(0);
		glBindTexture(GL_TEXTURE_2D, 0);

		
	}
}

namespace Graphene {


	/*struct Character
	{
		unsigned int TextureID;
		glm::ivec2 Size;
		glm::ivec2 Bearing;
		unsigned int Advance;
	};*/


	void framebufferSizeCallback(GLFWwindow* window, int width, int height);

	GrapheneColor::GrapheneColor(float red, float green, float blue)
	{
		this->r = red;
		this->g = green;
		this->b = blue;
	}
	GrapheneColor Grey(45.0f / 255.0f, 51.0f / 255.0f, 41.0f / 255.0f);
	GrapheneColor ambient_background_1(188.0f / 255.0f, 188.0f / 255.0f, 138.0f / 255.0f);



	

	Application::Application(int argc, char* argv[]) {
		if (!init(argc, argv))
		{
			std::cout << "Cannot initialize" << std::endl;
		}

		glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
		glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
		glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);



		//read logic here
		main = glfwCreateWindow(1280, 720, "App", NULL, NULL);
		if (!main)
		{
			std::cerr << "Could not make window" << std::endl;
			glfwTerminate();
			return;
		}
		glfwMakeContextCurrent(main);
		glfwSetFramebufferSizeCallback(main, framebufferSizeCallback);


		
		if (!gladLoadGL())
		{

			std::cerr << "Could not initialize GLAD" << std::endl;
			glfwTerminate();
			return;
		}


		int afbw, afbh;
		glfwGetFramebufferSize(main, &afbw, &afbh);
		glViewport(0, 0, afbw, afbh);
		SetVsync(true);
		
		//glEnable(GLFW_OPENGL_DEBUG_CONTEXT)
	};

	static void setClearColor(const GrapheneColor& color)
	{
		glClearColor(color.r, color.g, color.b, 1.0f);
	}



void framebufferSizeCallback(GLFWwindow* window, int width, int height)
{
	glViewport(0, 0, width, height);
}
void Application::onUpdate() {
	glfwSwapBuffers(main);
	glfwPollEvents();
}

void Application::GHClose() {
	glfwDestroyWindow(main);
	glfwTerminate();
}

GLFWwindow* Application::getWindow()
{
	return main;
}

void Application::GHStart() {
	std::cout << "Program Start" << std::endl;
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	int width, height;
	glfwGetWindowSize(main, &width, &height);

	//set uniform for projection matrix
	glm::mat4 projection = glm::ortho(0.0f, (float)width, 0.0f, (float)height);

	//grey.r, grey.g, grey.b
	//top rect
	GLfloat vertices[] =
	{
		//	Coordinates	//	Colors	//	Normals
		width, height, 0.0f,			Grey.r, Grey.g, Grey.b,
		width, height * 0.9f, 0.0f,		Grey.r, Grey.g, Grey.b,
		0.0f, height, 0.0f,				Grey.r, Grey.g, Grey.b,
		0.0f, height * 0.9f, 0.0f,		Grey.r, Grey.g, Grey.b,


		width, height * 0.1f, 0.0f,		Grey.r, Grey.g, Grey.b,
		width, 0.0f, 0.0f,				Grey.r, Grey.g, Grey.b,
		0.0f, height * 0.1f, 0.0f,		Grey.r, Grey.g, Grey.b,
		0.0f, 0.0f, 0.0f,				Grey.r, Grey.g, Grey.b,


	};
	GLuint indices[] =
	{
		0, 1, 2,
		1, 2, 3,
		4, 5, 6,
		5, 6, 7
	};
	VAO VAO1;
	VAO1.Bind();

	VBO VBO1(vertices, sizeof(vertices));
	EBO EBO1(indices, sizeof(indices));


	VAO1.LinkAttrib(VBO1, 0, 3, GL_FLOAT, 6 * sizeof(float), (void*)0);
	VAO1.LinkAttrib(VBO1, 1, 3, GL_FLOAT, 6 * sizeof(float), (void*)(3 * sizeof(float)));
	
	VAO1.Unbind();
	VBO1.Unbind();
	EBO1.Unbind();

	VAO VAOTXT;

	Shader textShader("../Graphene/src/glyphs.vert", "../Graphene/src/glyphs.frag");
	textShader.Activate();
	
	glUniformMatrix4fv(glGetUniformLocation(textShader.ID, "projection"), 1, GL_FALSE, glm::value_ptr(projection));
	
	setClearColor(ambient_background_1);

	
	std::vector<Glyph> glpyhs = loadFont("../Graphene/fonts/UltronRegular-ARwBM.ttf", 40);

	
	Shader objectShader("../Graphene/src/planar.vert", "../Graphene/src/planar.frag");
	objectShader.Activate();
	std::cout << glGetUniformLocation(objectShader.ID, "projection") << std::endl;
	glUniformMatrix4fv(glGetUniformLocation(objectShader.ID, "projection"), 1, GL_FALSE, glm::value_ptr(projection));

	
	//while (!glfwWindowShouldClose(main))
	//{

		glClear(GL_COLOR_BUFFER_BIT);

		VAO1.Bind();
		objectShader.Activate();
		
		glDrawElements(GL_TRIANGLES, sizeof(indices) / sizeof(int), GL_UNSIGNED_INT, 0);
		VAOTXT.Bind();
		textShader.Activate();
		//TODO(fix error here)
		renderText(textShader, "Graphene", VAOTXT, glpyhs, (float)width * 0.5, (float)height * 0.5, 1.0f, glm::vec3(1.0f, 1.0f, 1.0f));

			onUpdate();
		//}

		VAO1.Delete();   
		VBO1.Delete();
		EBO1.Delete();

	}

	bool Application::init(int argc = 0, char* argv[] = nullptr) {
		if (!glfwInit()) {
			std::cout << "INIT ERROR";
			return false;
		}
		return true;
	}
	void Application::SetVsync(bool set) {
		if (set) { glfwSwapInterval(1); }
	}


	//assumes plot
	void Application::GHSetupWindow(Window* Win_Title, uint16_t WinType, Dimensions* label_ptr, int PlotType) {
		int width = Win_Title->GetWidth();
		int height = Win_Title->GetHeight();

		//maybe just resize/restructure main window???
		//main = glfwCreateWindow(width, height, Win_Title->GetTitle(), NULL, NULL);

	}
	//assumes not plot
	void Application::GHSetupWindow(Window* Win_Title, uint16_t WinType, Dimensions* label_ptr) {
		int width = Win_Title->GetWidth();
		int height = Win_Title->GetHeight();

		main = glfwCreateWindow(width, height, Win_Title->GetTitle(), NULL, NULL);
	}

	void Application::GHSet2D(const char* PlotType, const char* X_Label, const char* Y_Label,
		const char* State,
		int ulx, int uly, SCALAR X_Scale, SCALAR Y_Scale,
		int X_Auto_Rescale, int Y_Auto_Rescale, SCALAR X_Min,
		SCALAR X_Max, SCALAR Y_Min, SCALAR Y_Max)
	{
		int             plottype;

		Dimensions* label_ptr = new Dimensions();

		//label_ptr = SetupLabelStruct();
		label_ptr->Y_Label = _strdup(Y_Label);
		label_ptr->Y_Min = Y_Min;
		label_ptr->Y_Max = Y_Max;
		label_ptr->Y_Scale = Y_Scale;
		label_ptr->Y_Auto_Rescale = Y_Auto_Rescale;

		label_ptr->X_Label = _strdup(X_Label);
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
		//GHSetupWindow(Y_Label, ulx, uly, WINDOW2D, label_ptr, plottype);
	}


	Module Application::GHModCreatePlot(const char* theInputFile) {

		//Search module directory for this file first

		FILE* inputfile = new FILE;
		FILE* tmp = new FILE;
		ArrayList<ArrayList<verts>> v;
		//convert plots to vertices?
		ArrayList<plot>* thePlots = new ArrayList<plot>;
		char line[512];
		char filename[512];
		char lastline[512];
		int xloc, yloc;
		//XInit(argc, argv, &theTime);
		if ((fopen_s(&inputfile, theInputFile, "r")) == NULL) {
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

		while (fgets(line, 511, inputfile)) {

			sscanf_s(line, "%s %d %d", filename, &xloc, &yloc);
			//	 sscanf(line,"%d %d",&xloc,&yloc);
			// 

			xloc = std::max(xloc, 1); //xloc = min(700, xloc);
			yloc = std::max(yloc, 1); //yloc = min(700, yloc);
			if (!strcmp(filename, "END")) break;
			std::vector<char> inV = file_to_buffer(tmp, 0, true);

			inputType t = determineFType(&tmp, filename, 0);
			switch (t)
			{
			case inputType::LINE_PLOT:
				/*LinePlot curPlot(inV, xloc, yloc);
				ArrayList<structureData> temp = curPlot.getStructures();
				for (int i = 0; i < temp.nItems; i++) {
					v.add(&plot_structure_to_vertices(temp[i], TWOD));
				}*/
				break;
			case inputType::SCATTER_PLOT:
				/*ScatterPlot curPlot(inV, xloc, yloc);
				ArrayList<structureData> temp = curPlot.getStructures();
				for (int i = 0; i < temp.nItems; i++) {
					v.add(&plot_structure_to_vertices(temp[i], TWOD));
				}*/
				break;
			case inputType::SURFACE_PLOT:
				/*SurfacePlot curPlot(inV, xloc, yloc);
				ArrayList<structureData> temp = curPlot.getStructures();
				for (int i = 0; i < temp.nItems; i++) {
					v.add(&plot_structure_to_vertices(temp[i], TWOD));
				}*/
				break;
			case inputType::VECTOR_PLOT:
				/*VectorPlot curPlot(inV, xloc, yloc);
				ArrayList<structureData> temp = curPlot.getStructures();
				for (int i = 0; i < temp.nItems; i++) {
					v.add(&plot_structure_to_vertices(temp[i], TWOD));
				}*/
				break;
			default:
				fprintf(stderr, "Unsupported plot type in file %s\n", filename);
			}
			fclose(tmp);
		}
		//create file using filename
		const char* t = *theInputFile + ".mod";
		Module m(t, v);
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
			//ColorCode_On(theNewWindow);
			//Grid_On(theNewWindow);
			//win->paint_function = Paint_ThreeD_Window;
			//win->print_function = PostScript_ThreeD_Window;
			//win->ascii_print_function = Ascii_ThreeD_Window;
			//win->xgrafix_print_function = Bin_ThreeD_Window;
			if (PlotType & 0b100) {}
			//	Set_X_Log(theNewWindow);
			if (PlotType & 0b010) {}
				//Set_Y_Log(theNewWindow);
			if (PlotType & 0b001) {}
				//Set_Z_Log(theNewWindow);
		}
		else if (WinType == WINDOW_VECTOR_PLOT) {
			//VectorWindow win(mainProps, PlotType);
			//win->paint_function = Paint_Vector_Window;
			//win->print_function = PostScript_Vector_Window;
			//win->ascii_print_function = Ascii_Vector_Window;
			//win->xgrafix_print_function = Bin_Vector_Window;
		}
		else {
			TWODWindow win(mainProps, PlotType);
			//theNewWindow->paint_function = Paint_Window;
			//theNewWindow->print_function = PostScriptOpenWindow;
			//theNewWindow->ascii_print_function = Ascii_TwoD_Window;
			//theNewWindow->xgrafix_print_function = Bin_TwoD_Window;
			if (PlotType & 0b10) {}
				//Set_X_Log(theNewWindow);
				if (PlotType & 0b01) {}
				//Set_Y_Log(theNewWindow);
		}

		if (label_ptr->Z_Auto_Rescale) {}
			//Z_AutoRescale(theNewWindow);
			if (label_ptr->Y_Auto_Rescale) {}
			//Y_AutoRescale(theNewWindow);
			if (label_ptr->X_Auto_Rescale) {}
			//X_AutoRescale(theNewWindow);

		return;

	}


	
};