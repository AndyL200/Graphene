#pragma once

#ifdef GH_PLATFORM_WINDOWS


extern Graphene::Application* CreateApplication();


	int main() {
		auto app = CreateApplication();
		app->GHStart();

		app->GHClose();
		return 0;
	}

#endif