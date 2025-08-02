#pragma once
#include "Application.h"
#include "Graphene.h"

extern Graphene::Application* CreateApplication();


	int main() {
		Graphene::Application* app = CreateApplication();
		app->GHStart();

		app->GHClose();
		return 0;
	}

