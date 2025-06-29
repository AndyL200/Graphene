#ifndef GRAPENE_PLOTS
#define GRAPHENE_PLOTS

#include "Graphene.h"
#include "Lists.h"
#include "Labels.h"

static int MAX_LIN = 12000;

struct structureData {
	int n;
	int fillFlag;
	int lineColor;
	int fillColor;
	SCALAR* points;
};


enum PlotCategory { // compiling xgmovie was throwing error for unknown LINE_PLOT; JK 2019-01-14;
	None = 0,
	LinePlotCategory = BIT(0),				
	ScatterPlotCategory = BIT(1),
	VectorPlotCategory = BIT(2),
	SurfacePlotCategory = BIT(3),
	Scatter3DCategory = BIT(4),
	Irregular3DCategory = BIT(5)
};


class plot 
{
protected:
	std::vector<uint8_t>& buffer;
	PlotDimensions dim;
	ArrayList<structureData> structures;
public:
	int xloc;
	int yloc;

	virtual void updatePlot(double time) = 0;  // bring the plot up to this time
	plot(std::vector<uint8_t>& buffer, int xloc, int yloc);
	virtual ~plot() {};
};

struct One_D_plot_data {
public:
	SCALAR* x;
	SCALAR* y;
	int n;
	int color;
	SCALAR time;
	One_D_plot_data() {};
};


class One_D_plot : public plot {

protected:
	ArrayList< ArrayList<One_D_plot_data> > graphdata;
	ArrayList<One_D_plot_data>* current;
	//ArrayList< One_D_plot_data> current_data;
public:
	virtual void updatePlot(double time);
	One_D_plot(VectorStream& in, int xloc, int yloc);
};

class ScatterPlot : public One_D_plot {
public:
	ScatterPlot(VectorStream& in, int xloc, int yloc);
};

class LinePlot : public One_D_plot {
public:
	LinePlot(VectorStream& in, int xloc, int yloc);
};



class SurfacePlotData {
public:
	SCALAR** z;
	SCALAR time;
};

class SurfacePlot : public plot {
public:
	int n, m;
	SCALAR* x;
	SCALAR* y;
	SCALAR** z;
	ArrayList< SurfacePlotData > graphdata;
	ListIter<SurfacePlotData>* current;
	virtual void updatePlot(double time);
	SurfacePlot(char* _filename, int xloc, int yloc);
};

class VectorPlotData {
public:
	SCALAR** z;
	SCALAR** w;
	SCALAR time;
};

class VectorPlot : public plot {
public:
	int n, m;
	SCALAR* x;
	SCALAR* y;
	SCALAR** z;
	SCALAR** w;
	ArrayList< VectorPlotData > graphdata;
	ListIter<VectorPlotData>* current;
	virtual void updatePlot(double time);
	VectorPlot(char* _filename, int xloc, int yloc);
};


#endif