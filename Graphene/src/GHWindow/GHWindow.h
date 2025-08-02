#ifndef GRAPHENE_WINDOW
#define GRAPHENE_WINDOW


#include "Core.h"
#include "ShaderHead.h"


namespace Graphene {

	enum WindowType : uint16_t {
		WINDOW2D = BIT(0),
		WINDOW3D = BIT(1),
		WINDOW_PLOT = BIT(2),
		WINDOW_ONE_D_PLOT = BIT(3),
		WINDOW_SCATTER_PLOT = BIT(4),
		WINDOW_VECTOR_PLOT = BIT(5),
		WINDOW_FREE = BIT(6)
	};

	
	struct WindowProps
	{
		std::string Title;
		unsigned int Width;
		unsigned int Height;
		WindowProps operator=(const WindowProps& w) {
			this->Title = w.Title;
			this->Width = w.Width;
			this->Height = w.Height;
			return *this;
		}
		WindowProps(const std::string& title = "Graphene Graphics", unsigned int width = 1280, unsigned int height = 720) : Title(title), Width(width), Height(height)
		{

		}

	};

	

	class Window {
	public:
		//using EventCallbackFn = std::function<void(Event&)>;

		Window(WindowProps& props) {
			this->props = props;
		}

		virtual ~Window() {}

		void plotInit() {}

		virtual unsigned int GetWidth() = 0;
		virtual unsigned int GetHeight() = 0;
		virtual const char* GetTitle() = 0;
		virtual void OnUpdate() = 0;

		static GLFWwindow* Create() {};

	protected:
		WindowProps props;
	};

	class TWODWindow : public Window {
	public:
		TWODWindow(WindowProps props, int type);

		GLFWwindow* Create();

		//theory, inline may help keep the context defined in Application
		void inline plotInit(GLFWwindow* window);

		void freeInit(GLFWwindow* window);
		void SetVertices(std::vector<float>& v);

		GLfloat* GetVertices();

		unsigned int GetWidth() override;
		unsigned int GetHeight() override;
		const char* GetTitle() override;

		void OnUpdate() override;

		int type;
		GLfloat* vertices;
	};



	class THREEDWindow : public Window {
	public:
		THREEDWindow(WindowProps props, int type);
		unsigned int GetWidth() override;
		unsigned int GetHeight() override;
		const char* GetTitle() override;
		GLFWwindow* Create();

		void plotInit(GLFWwindow* window);
		

		void freeInit(GLFWwindow* window);

		void OnUpdate() override;

		int type;
		GLfloat* vertices;
	};

	/*static GLFWwindow* CreateWindow(const WindowProps& props = WindowProps(), const int type = WINDOW_FREE)
	{
		if (type & WINDOW2D) {
			TWODWindow window(props, type);
			return window.Create();
		}
		else if (type & WINDOW3D) {
			THREEDWindow window(props, type);
			return window.Create();
		}
		else {
			return nullptr;
		}
	}*/

}




#endif