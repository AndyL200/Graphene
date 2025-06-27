#ifndef GRAPENE_PLOTS
#define GRAPHENE_PLOTS

#include "Graphene.h"

struct structureData {
	int n;
	int fillFlag;
	int lineColor;
	int fillColor;
	SCALAR* points;
};

class plot 
{
protected:
	std::string filename;
	LabelType label;
	List<StructureData> structures;
public:
	int xloc;
	int yloc;

	virtual void updatePlot(double time) = 0;  // bring the plot up to this time
	plot(std::string filename, int xloc, int yloc);
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
	List< List<One_D_plot_data> > graphdata;
	ListIter<List<One_D_plot_data> >* current;
	List< One_D_plot_data> current_data;
public:
	virtual void updatePlot(double time);
	One_D_plot(char* _filename, int xloc, int yloc);
};

class ScatterPlot : public One_D_plot {
public:
	ScatterPlot(char* _filename, int xloc, int yloc);
};

class LinePlot : public One_D_plot {
public:
	LinePlot(char* _filename, int xloc, int yloc);
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
	List< SurfacePlotData > graphdata;
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
	List< VectorPlotData > graphdata;
	ListIter<VectorPlotData>* current;
	virtual void updatePlot(double time);
	VectorPlot(char* _filename, int xloc, int yloc);
};


#endif