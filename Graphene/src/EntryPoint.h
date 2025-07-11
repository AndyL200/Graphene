#pragma once
#include "Application.h"
#include "Graphene.h"


namespace Graphene {
	int main() {
		Application app;
		app.GHStart();

		app.GHClose();
		return 0;
	}

}