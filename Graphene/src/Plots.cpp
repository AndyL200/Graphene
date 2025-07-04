#include "Plots.h"
#include <xutility>

plot::plot(std::vector<uint8_t>& buffer, int _xloc, int _yloc) {
	
	xloc = _xloc;
	yloc = _yloc;

	fclose(fp);
}


PlotDimensions ReadLabel(VectorStream& in) {
	int type;
	char buf[100];
	char* labels;
	int l;
	PlotDimensions* label = (PlotDimensions*)malloc(sizeof(PlotDimensions));
	SCALAR scale, min, max;
	//Read 1st string
	l = GHRead<int>(in);  //size
	labels = (char*)malloc((l + 1) * sizeof(char));
	if (l > 0) {
		in.read(labels, l);
		labels[l] = '\0';
		//this could throw off the vector stream pointer
		//make a read specifically for pointer values?
		label->X_Label = labels;
	}

	//Read 2st string
	l = GHRead<int>(in);  //size
	labels = (char*)malloc((l + 1) * sizeof(char));
	if (l > 0) {
		in.read(labels, l);
		labels[l] = '\0';
		label->Y_Label = labels;
	}
	//Read 3st string
	l = GHRead<int>(in);  //size
	labels = (char*)malloc((l + 1) * sizeof(char));
	if (l > 0) {
		in.read(labels, l);
		labels[l] = '\0';
		label->Z_Label = labels;
	}
	// read the info on the X axis
	scale = GHRead<SCALAR>(in);
	label->X_Scale = 1;//scale;
	min = GHRead<SCALAR>(in);
	label->X_Min = min;
	max = GHRead<SCALAR>(in);
	label->X_Max = max;
	label->X_Auto_Rescale = GHRead<SCALAR>(in);


	// read the info on the Y axis
	scale = GHRead<SCALAR>(in);
	label->Y_Scale = 1;//scale;
	min = GHRead<SCALAR>(in);
	label->Y_Min = min;
	max = GHRead<SCALAR>(in);
	label->Y_Max = max;
	label->Y_Auto_Rescale = GHRead<SCALAR>(in);


	// read the info on the Z axis
	scale = GHRead<SCALAR>(in);
	label->Z_Scale = 1;//scale;
	min = GHRead<SCALAR>(in);
	label->Z_Min = min;
	max = GHRead<SCALAR>(in);
	label->Z_Max = max;
	label->Z_Auto_Rescale = GHRead<SCALAR>(in);

	return *label;
}


One_D_plot::One_D_plot(VectorStream& in, int xloc, int yloc) : plot(buffer, xloc, yloc) {
	int type;
	ArrayList<ArrayList<One_D_plot_data> > tmp;
	const int count = GHRead<int>(in);
	in.clear();
	in.seekg(0, std::ios::beg);
	//  for(i=1;openFile(&fp,filename,i)!=-1;i++) {
	for (int i = 0; i < count; i++) {
		ArrayList<One_D_plot_data>* in_data = new ArrayList<One_D_plot_data>;
		SCALAR the_time;
		type = GHRead<int>(in);
		dim = ReadLabel(in);
		the_time = GHRead<SCALAR>(in);
		//issue ArrayList cannot be read
		do {
			One_D_plot_data* dat = new One_D_plot_data();
			dat->n = GHRead<int>(in);
			dat->color = GHRead<int>(in);
			int Skip = 1;
			
			if (dat->n == -1) {  // the exit condition 
				tmp.add(in_data);  //add this list of data to the graph list
				break;  // go to the next file
			}

			if (dat->n > MAX_LIN) {
				// we need to contract the dataset to fit into our
				// memory buffer.  We also must read the entire dataset.
				// what we do is stick dat->n/MAX_LIN items directly
				// into the same array location, the last one wins.
				Skip = dat->n / MAX_LIN + 1;

				if (Skip < 2) Skip = 2;
				//		  printf("\n dat->n %d, Skip %d ",dat->n, Skip);
					printf("Warning, some datasets are too large.  Threading them. (%d)\n", Skip);
			}

			dat->x = (SCALAR*)calloc(max(1, (dat->n + Skip) / Skip), sizeof(SCALAR));
			dat->y = (SCALAR*)calloc(max(1, (dat->n + Skip) / Skip), sizeof(SCALAR));
			int idx = 0;
			for (int j = 0; j < dat->n; j++) {
				if (j % Skip != 0)
					continue;
				SCALAR xval = GHRead<SCALAR>(in);
				SCALAR yval = GHRead<SCALAR>(in);
				dat->x[idx] = xval;
				dat->y[idx] = yval;
				++idx;

			}
			dat->n = idx;
			in_data->add(dat);
		} while (in.eof());  // actually a dummy, reading should stop earlier
	}
	// reverse the list
	graphdata = tmp;

}

void One_D_plot::updatePlot(double time) {
	//if (time == minTime) current->restart();

	ArrayList<One_D_plot_data>* thisdata = current;
	while (!current->Done() && current->current()->head->data->time <= time)
	{
		thisdata = current->current();
		(*current)++;
	}
	if (thisdata) {
		ListIter<One_D_plot_data> walk(*thisdata);
		ListIter<One_D_plot_data> walk2(current_data);
		for (walk.restart(), walk2.restart(); !walk.Done(); walk++, walk2++) {
			One_D_plot_data* thisgraph = walk.current();
			One_D_plot_data* destgraph = walk2.current();
			destgraph->n = thisgraph->n;
			memcpy(destgraph->x, thisgraph->x, destgraph->n * sizeof(SCALAR));
			memcpy(destgraph->y, thisgraph->y, destgraph->n * sizeof(SCALAR));
		}
	}
}




ScatterPlot::ScatterPlot(char* filename, int xloc, int yloc) :One_D_plot(filename, xloc, yloc) {
	XGSet2D("linlin", dim->X_Label, dim->Y_Label, "open", xloc, yloc, dim->X_Scale,
		dim->Y_Scale, dim->X_Auto_Rescale, dim->Y_Auto_Rescale,
		dim->X_Min, dim->X_Max, dim->Y_Min, dim->Y_Max);
	ListIter<One_D_plot_data> walk(current_data);
	for (walk.restart(); !walk.Done(); walk++)
		XGScat2D(walk.current()->x, walk.current()->y, &(walk.current()->n), walk.current()->color);
	ListIter<StructureData> walk2(structures);
	for (walk2.restart(); !walk2.Done(); walk2++)
		XGStructureArray(walk2.current()->n, (STRUCT_FILL)walk2.current()->fillFlag,
			walk2.current()->lineColor,
			walk2.current()->fillColor, walk2.current()->points);
}

LinePlot::LinePlot(char* filename, int xloc, int yloc) :One_D_plot(filename, xloc, yloc) {
	XGSet2D("linlin", dim->X_Label, dim->Y_Label, "open", xloc, yloc, dim->X_Scale,
		dim->Y_Scale, dim->X_Auto_Rescale, dim->Y_Auto_Rescale,
		dim->X_Min, dim->X_Max, dim->Y_Min, dim->Y_Max);
	ListIter<One_D_plot_data> walk(current_data);
	for (walk.restart(); !walk.Done(); walk++)
		XGCurve(walk.current()->x, walk.current()->y, &(walk.current()->n), walk.current()->color);
	ListIter<StructureData> walk2(structures);
	for (walk2.restart(); !walk2.Done(); walk2++)
		XGStructureArray(walk2.current()->n, (STRUCT_FILL)walk2.current()->fillFlag,
			walk2.current()->lineColor,
			walk2.current()->fillColor, walk2.current()->points);
}

SurfacePlot::SurfacePlot(char* filename, int xloc, int yloc) : plot(filename, xloc, yloc) {
	FILE* fp;
	SCALAR the_time;
	List<SurfacePlotData> tmp;
	int count, type;

	x = 0; y = 0;
	//  for(int i=1;openFile(&fp,filename,i)!=-1;i++) {
	openFile(&fp, filename, 0);
	count = getCount(filename);
	rewind(fp);
	for (int i = 0; i < count; i++) {
		SurfacePlotData* thisgraph = new SurfacePlotData;
		XGRead(&type, sizeof(int), 1, fp, "int");
		ReadLabel(fp);
		XGRead(&(thisgraph->time), sizeof(SCALAR), 1, fp, SCALAR_CHAR);
		maxTime = max(maxTime, thisgraph->time);
		minTime = min(minTime, thisgraph->time);
		XGRead(&n, sizeof(int), 1, fp, "int");
		XGRead(&m, sizeof(int), 1, fp, "int");
		if (!x) {
			x = (SCALAR*)malloc(m * sizeof(SCALAR));
			y = (SCALAR*)malloc(n * sizeof(SCALAR));
			z = (SCALAR**)malloc(m * sizeof(SCALAR*));
			for (int j = 0; j < m; j++) {
				z[j] = (SCALAR*)calloc(n, sizeof(SCALAR));
				//	  memset(z[j],0,n * sizeof(SCALAR));
			}

		}
		// read the x and y axes
		XGRead(x, sizeof(SCALAR), m, fp, SCALAR_CHAR);
		XGRead(y, sizeof(SCALAR), n, fp, SCALAR_CHAR);


		thisgraph->z = (SCALAR**)malloc(m * sizeof(SCALAR*));
		for (int j = 0; j < m; j++) {
			thisgraph->z[j] = (SCALAR*)malloc(n * sizeof(SCALAR));
			XGRead(thisgraph->z[j], sizeof(SCALAR), n, fp, SCALAR_CHAR);
		}
		tmp.add(thisgraph);
	}
	fclose(fp);

	graphdata = tmp;  //this will reverse tmp
	current = new ListIter<SurfacePlotData>(graphdata);
	XGSet3D("linlinlin", dim->X_Label, dim->Y_Label, dim->Z_Label, 45.0, 225.0,
		"open", xloc, yloc, dim->X_Scale, dim->Y_Scale, dim->Z_Scale,
		dim->X_Auto_Rescale, dim->Y_Auto_Rescale, dim->Z_Auto_Rescale,
		dim->X_Min, dim->X_Max, dim->Y_Min, dim->Y_Max, dim->Z_Min, dim->Z_Max);
	XGSurf(x, y, z, &m, &n, 1);

}

void SurfacePlot::updatePlot(double time) {

	if (time == minTime) current->restart();

	SurfacePlotData* thisdata = current->current();
	while (!current->Done() && current->current()->time <= time)
	{
		thisdata = current->current();
		(*current)++;
	}
	if (thisdata)
		for (int j = 0; j < m; j++)
			memcpy(z[j], thisdata->z[j], n * sizeof(SCALAR));
}



VectorPlot::VectorPlot(char* filename, int xloc, int yloc) : plot(filename, xloc, yloc) {
	FILE* fp;
	SCALAR the_time;
	List<VectorPlotData> tmp;
	int count, type;
	x = 0; y = 0;

	//  for(int i=1;openFile(&fp,filename,i)!=-1;i++) {
	openFile(&fp, filename, 0);
	count = getCount(filename);
	rewind(fp);
	for (int i = 0; i < count; i++) {
		VectorPlotData* thisgraph = new VectorPlotData;
		XGRead(&type, sizeof(int), 1, fp, "int");
		ReadLabel(fp);
		XGRead(&(thisgraph->time), sizeof(SCALAR), 1, fp, SCALAR_CHAR);
		maxTime = max(maxTime, thisgraph->time);
		minTime = min(minTime, thisgraph->time);
		XGRead(&n, sizeof(int), 1, fp, "int");
		XGRead(&m, sizeof(int), 1, fp, "int");
		if (!x) {
			x = (SCALAR*)malloc(m * sizeof(SCALAR));
			y = (SCALAR*)malloc(n * sizeof(SCALAR));
			z = (SCALAR**)malloc(m * sizeof(SCALAR*));
			w = (SCALAR**)malloc(m * sizeof(SCALAR*));
			for (int j = 0; j < m; j++) {
				z[j] = (SCALAR*)calloc(n, sizeof(SCALAR));
				w[j] = (SCALAR*)calloc(n, sizeof(SCALAR));
			}

		}
		// read the x and y axes
		XGRead(x, sizeof(SCALAR), m, fp, SCALAR_CHAR);
		XGRead(y, sizeof(SCALAR), n, fp, SCALAR_CHAR);


		thisgraph->z = (SCALAR**)malloc(m * sizeof(SCALAR*));
		thisgraph->w = (SCALAR**)malloc(m * sizeof(SCALAR*));
		for (int j = 0; j < m; j++) {
			thisgraph->w[j] = (SCALAR*)malloc(n * sizeof(SCALAR));
			thisgraph->z[j] = (SCALAR*)malloc(n * sizeof(SCALAR));
			XGRead(thisgraph->w[j], sizeof(SCALAR), n, fp, SCALAR_CHAR);
			XGRead(thisgraph->z[j], sizeof(SCALAR), n, fp, SCALAR_CHAR);
		}
		tmp.add(thisgraph);
	}
	fclose(fp);

	graphdata = tmp;  //this will reverse tmp
	current = new ListIter<VectorPlotData>(graphdata);
	XGSetVec("vecvec", label->X_Label, label->Y_Label, label->Z_Label,
		"open", xloc, yloc, label->X_Scale, label->Y_Scale,
		label->X_Auto_Rescale, label->Y_Auto_Rescale,
		label->X_Min, label->X_Max, label->Y_Min, label->Y_Max);
	XGVector(x, y, w, z, &m, &n, 1);
	ListIter<StructureData> walk2(structures);
	for (walk2.restart(); !walk2.Done(); walk2++)
		XGStructureArray(walk2.current()->n, (STRUCT_FILL)walk2.current()->fillFlag,
			walk2.current()->lineColor,
			walk2.current()->fillColor, walk2.current()->points);
}

void VectorPlot::updatePlot(double time) {

	if (time == minTime) current->restart();

	VectorPlotData* thisdata = current->current();
	while (!current->Done() && current->current()->time <= time)
	{
		thisdata = current->current();
		(*current)++;
	}
	if (thisdata)
		for (int j = 0; j < m; j++) {
			memcpy(z[j], thisdata->z[j], n * sizeof(SCALAR));
			memcpy(w[j], thisdata->w[j], n * sizeof(SCALAR));
		}
}

List<plot>* thePlots;

int main(int argc, char** argv) {
	FILE* inputfile;
	FILE* tmp;
	thePlots = new List<plot>;
	char line[512];
	char filename[512];
	char lastline[512];
	int xloc, yloc;
	XGInit(argc, argv, &theTime);
	if ((inputfile = fopen(theInputFile, "r")) == NULL) {
		fprintf(stderr, "Cannot open inputfile.  Exiting.");
		exit(1);
	}
	theTime = 0;
	maxTime = 0;
	minTime = 1e9;
	fgets(line, 511, inputfile);
#ifdef XG_SCALAR_DOUBLE
	sscanf(line, "%lf %d", &thetimestep, &MAX_LIN);
#else
	sscanf(line, "%f %d", &thetimestep, &MAX_LIN);
#endif
	if (MAX_LIN <= 0) MAX_LIN = 12000;

	while (!feof(inputfile)) {
		fgets(line, 511, inputfile);
		sscanf(line, "%s %d %d", filename, &xloc, &yloc);
		//	 sscanf(line,"%d %d",&xloc,&yloc);
		xloc = max(xloc, 1); xloc = min(700, xloc);
		yloc = max(yloc, 1); yloc = min(700, yloc);
		if (!strcmp(filename, "END")) break;
		switch (openFile(&tmp, filename, 0))
		{
		case LINE_PLOT:
			thePlots->add(new LinePlot(filename, xloc, yloc));
			break;
		case SCATTER_PLOT:
			thePlots->add(new ScatterPlot(filename, xloc, yloc));
			break;
		case SURFACE_PLOT:
			thePlots->add(new SurfacePlot(filename, xloc, yloc));
			break;
		case VECTOR_PLOT:
			thePlots->add(new VectorPlot(filename, xloc, yloc));
			break;
		default:
			fprintf(stderr, "Unsupported plot type in file %s\n", filename);
		}
		fclose(tmp);
	}
	theTime = minTime;

	ListIter<plot> walk(*thePlots);
	for (walk.restart(); !walk.Done(); walk++) {
		walk.current()->updatePlot(theTime);
	}
	printf("\nStart time: %g  End time: %g\n", minTime, maxTime);
	XGStart();
	return 0;
}

void XGMainLoop() {

	//  cause the graphs to cycle
	if (theTime > maxTime) theTime = minTime;

	ListIter<plot> walk(*thePlots);
	for (walk.restart(); !walk.Done(); walk++) {
		walk.current()->updatePlot(theTime);
	}
	theTime += thetimestep;

}

void Dump(const char*) {}
void Quit() {}
