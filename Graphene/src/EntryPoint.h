#pragma once
#include "Graphene.h"
#include "Application.h"

namespace Graphene {
	int main() {
		Application app;
		app.GHStart();

		app.GHClose();
		return 0;
	}

}