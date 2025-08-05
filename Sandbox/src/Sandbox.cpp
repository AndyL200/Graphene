#include "Graphene.h"


class Sandbox : public Graphene::Application
{
public:
	Sandbox() 
	{
	
	}
	~Sandbox()
	{

	}
};

	Graphene::Application* CreateApplication()
	{
		return new Sandbox();
	}


