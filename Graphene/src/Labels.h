#ifndef GRAPHENE_LABELS
#define GRAPHENE_LABELS

struct Dimensions {
	double          Z_Scale = 0;
	double          Z_Min = 0;
	double          Z_Max = 0;
	char* Z_Label = nullptr;
	int             Z_Auto_Rescale = 0;
	double          Y_Scale = 0;
	double          Y_Offset = 0;
	double          Y_Min = 0;
	double          Y_Max = 0;
	char* Y_Label = nullptr;
	int             Y_Auto_Rescale = 0;
	double          X_Scale = 0;
	double          X_Offset = 0;
	double          X_Min = 0;
	double          X_Max = 0;
	char* X_Label = nullptr;
	int             X_Auto_Rescale = 0;

	Dimensions() {}
};

#endif