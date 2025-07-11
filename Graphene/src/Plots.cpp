#include "Plots.h"



Dimensions ReadLabel(VectorStream& in) {
	int type;
	char buf[100];
	char* labels;
	int l;
	Dimensions* label = (Dimensions*)malloc(sizeof(Dimensions));
	SCALAR scale, min, max;
	//Read 1st string
	l = *GHRead<int>(in);  //size
	labels = (char*)malloc((l + 1) * sizeof(char));
	if (l > 0) {
		in.read(labels, l);
		labels[l] = '\0';
		//this could throw off the vector stream pointer
		//make a read specifically for pointer values?
		label->X_Label = labels;
	}

	//Read 2st string
	l = *GHRead<int>(in);  //size
	labels = (char*)malloc((l + 1) * sizeof(char));
	if (l > 0) {
		in.read(labels, l);
		labels[l] = '\0';
		label->Y_Label = labels;
	}
	//Read 3st string
	l = *GHRead<int>(in);  //size
	labels = (char*)malloc((l + 1) * sizeof(char));
	if (l > 0) {
		in.read(labels, l);
		labels[l] = '\0';
		label->Z_Label = labels;
	}
	// read the info on the X axis
	scale = *GHRead<SCALAR>(in);
	label->X_Scale = 1;//scale;
	min = *GHRead<SCALAR>(in);
	label->X_Min = min;
	max = *GHRead<SCALAR>(in);
	label->X_Max = max;
	label->X_Auto_Rescale = *GHRead<SCALAR>(in);


	// read the info on the Y axis
	scale = *GHRead<SCALAR>(in);
	label->Y_Scale = 1;//scale;
	min = *GHRead<SCALAR>(in);
	label->Y_Min = min;
	max = *GHRead<SCALAR>(in);
	label->Y_Max = max;
	label->Y_Auto_Rescale = *GHRead<SCALAR>(in);


	// read the info on the Z axis
	scale = *GHRead<SCALAR>(in);
	label->Z_Scale = 1;//scale;
	min = *GHRead<SCALAR>(in);
	label->Z_Min = min;
	max = *GHRead<SCALAR>(in);
	label->Z_Max = max;
	label->Z_Auto_Rescale = *GHRead<SCALAR>(in);

	return *label;
}


One_D_plot::One_D_plot(std::vector<char>& v, int xloc, int yloc) : plot(v, xloc, yloc) {
	int type;
	ArrayList<ArrayList<One_D_plot_data> > tmp;
	int count = *GHRead<int>(buffer);
	buffer.clear();
	buffer.seekg(0, std::ios::beg);
	//  for(i=1;openFile(&fp,filename,i)!=-1;i++) {
	for (int i = 0; i < count; i++) {
		ArrayList<One_D_plot_data>* in_data = new ArrayList<One_D_plot_data>;
		SCALAR the_time;
		type = *GHRead<int>(buffer);
		this->dim = ReadLabel(buffer);
		the_time = *GHRead<SCALAR>(buffer);
		//issue ArrayList cannot be read
		while(!buffer.eof()) {
			One_D_plot_data* dat = new One_D_plot_data();
			dat->n = *GHRead<int>(buffer);
			dat->color = *GHRead<int>(buffer);
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

			dat->x = (SCALAR*)calloc(std::max<size_t>(static_cast<size_t>(1), static_cast<size_t>((dat->n + Skip) / Skip)), sizeof(SCALAR));
			dat->y = (SCALAR*)calloc(std::max<size_t>(static_cast<size_t>(1), static_cast<size_t>((dat->n + Skip) / Skip)), sizeof(SCALAR));
			int idx = 0;
			for (int j = 0; j < dat->n; j++) {
				if (j % Skip != 0)
					continue;
				SCALAR xval = *GHRead<SCALAR>(buffer);
				SCALAR yval = *GHRead<SCALAR>(buffer);
				dat->x[idx] = xval;
				dat->y[idx] = yval;
				++idx;

			}
			dat->n = idx;
			in_data->add(dat);
		}  // actually a dummy, reading should stop earlier
	}
	// reverse the list
	graphdata = tmp;

}

//void One_D_plot::updatePlot(double time) {
//	//if (time == minTime) current->restart();
//
//	ArrayList<One_D_plot_data>* thisdata = current;
//	while (!current->Done() && current->current()->head->data->time <= time)
//	{
//		thisdata = current->current();
//		(*current)++;
//	}
//	if (thisdata) {
//		ListIter<One_D_plot_data> walk(*thisdata);
//		ListIter<One_D_plot_data> walk2(current_data);
//		for (walk.restart(), walk2.restart(); !walk.Done(); walk++, walk2++) {
//			One_D_plot_data* thisgraph = walk.current();
//			One_D_plot_data* destgraph = walk2.current();
//			destgraph->n = thisgraph->n;
//			memcpy(destgraph->x, thisgraph->x, destgraph->n * sizeof(SCALAR));
//			memcpy(destgraph->y, thisgraph->y, destgraph->n * sizeof(SCALAR));
//		}
//	}
//}




//ScatterPlot::ScatterPlot(vector<char>& v, int xloc, int yloc) :One_D_plot(v, xloc, yloc) {
//	XGSet2D("linlin", dim->X_Label, dim->Y_Label, "open", xloc, yloc, dim->X_Scale,
//		dim->Y_Scale, dim->X_Auto_Rescale, dim->Y_Auto_Rescale,
//		dim->X_Min, dim->X_Max, dim->Y_Min, dim->Y_Max);
//	ListIter<One_D_plot_data> walk(current_data);
//	for (walk.restart(); !walk.Done(); walk++)
//		XGScat2D(walk.current()->x, walk.current()->y, &(walk.current()->n), walk.current()->color);
//	ListIter<StructureData> walk2(structures);
//	for (walk2.restart(); !walk2.Done(); walk2++)
//		XGStructureArray(walk2.current()->n, (STRUCT_FILL)walk2.current()->fillFlag,
//			walk2.current()->lineColor,
//			walk2.current()->fillColor, walk2.current()->points);
//}

//LinePlot::LinePlot(vector<char>& v, int xloc, int yloc) :One_D_plot(v, xloc, yloc) {
//	XGSet2D("linlin", dim->X_Label, dim->Y_Label, "open", xloc, yloc, dim->X_Scale,
//		dim->Y_Scale, dim->X_Auto_Rescale, dim->Y_Auto_Rescale,
//		dim->X_Min, dim->X_Max, dim->Y_Min, dim->Y_Max);
//	ListIter<One_D_plot_data> walk(current_data);
//	for (walk.restart(); !walk.Done(); walk++)
//		XGCurve(walk.current()->x, walk.current()->y, &(walk.current()->n), walk.current()->color);
//	ListIter<StructureData> walk2(structures);
//	for (walk2.restart(); !walk2.Done(); walk2++)
//		XGStructureArray(walk2.current()->n, (STRUCT_FILL)walk2.current()->fillFlag,
//			walk2.current()->lineColor,
//			walk2.current()->fillColor, walk2.current()->points);
//}

//SurfacePlot::SurfacePlot(vector<char>& v, int xloc, int yloc) : plot(v, xloc, yloc) {
//	FILE* fp;
//	SCALAR the_time;
//	List<SurfacePlotData> tmp;
//	int count, type;
//
//	x = 0; y = 0;
//	//  for(int i=1;openFile(&fp,filename,i)!=-1;i++) {
//	openFile(&fp, filename, 0);
//	count = getCount(filename);
//	rewind(fp);
//	for (int i = 0; i < count; i++) {
//		SurfacePlotData* thisgraph = new SurfacePlotData;
//		XGRead(&type, sizeof(int), 1, fp, "int");
//		ReadLabel(fp);
//		XGRead(&(thisgraph->time), sizeof(SCALAR), 1, fp, SCALAR_CHAR);
//		maxTime = max(maxTime, thisgraph->time);
//		minTime = min(minTime, thisgraph->time);
//		XGRead(&n, sizeof(int), 1, fp, "int");
//		XGRead(&m, sizeof(int), 1, fp, "int");
//		if (!x) {
//			x = (SCALAR*)malloc(m * sizeof(SCALAR));
//			y = (SCALAR*)malloc(n * sizeof(SCALAR));
//			z = (SCALAR**)malloc(m * sizeof(SCALAR*));
//			for (int j = 0; j < m; j++) {
//				z[j] = (SCALAR*)calloc(n, sizeof(SCALAR));
//				//	  memset(z[j],0,n * sizeof(SCALAR));
//			}
//
//		}
//		// read the x and y axes
//		XGRead(x, sizeof(SCALAR), m, fp, SCALAR_CHAR);
//		XGRead(y, sizeof(SCALAR), n, fp, SCALAR_CHAR);
//
//
//		thisgraph->z = (SCALAR**)malloc(m * sizeof(SCALAR*));
//		for (int j = 0; j < m; j++) {
//			thisgraph->z[j] = (SCALAR*)malloc(n * sizeof(SCALAR));
//			XGRead(thisgraph->z[j], sizeof(SCALAR), n, fp, SCALAR_CHAR);
//		}
//		tmp.add(thisgraph);
//	}
//	fclose(fp);
//
//	graphdata = tmp;  //this will reverse tmp
//	current = new ListIter<SurfacePlotData>(graphdata);
//	XGSet3D("linlinlin", dim->X_Label, dim->Y_Label, dim->Z_Label, 45.0, 225.0,
//		"open", xloc, yloc, dim->X_Scale, dim->Y_Scale, dim->Z_Scale,
//		dim->X_Auto_Rescale, dim->Y_Auto_Rescale, dim->Z_Auto_Rescale,
//		dim->X_Min, dim->X_Max, dim->Y_Min, dim->Y_Max, dim->Z_Min, dim->Z_Max);
//	XGSurf(x, y, z, &m, &n, 1);
//
//}

//void SurfacePlot::updatePlot(double time) {
//
//	if (time == minTime) current->restart();
//
//	SurfacePlotData* thisdata = current->current();
//	while (!current->Done() && current->current()->time <= time)
//	{
//		thisdata = current->current();
//		(*current)++;
//	}
//	if (thisdata)
//		for (int j = 0; j < m; j++)
//			memcpy(z[j], thisdata->z[j], n * sizeof(SCALAR));
//}
//


//VectorPlot::VectorPlot(vector<char>& v, int xloc, int yloc) : plot(v, xloc, yloc) {
//	FILE* fp;
//	SCALAR the_time;
//	List<VectorPlotData> tmp;
//	int count, type;
//	x = 0; y = 0;
//
//	//  for(int i=1;openFile(&fp,filename,i)!=-1;i++) {
//	openFile(&fp, filename, 0);
//	count = getCount(filename);
//	rewind(fp);
//	for (int i = 0; i < count; i++) {
//		VectorPlotData* thisgraph = new VectorPlotData;
//		XGRead(&type, sizeof(int), 1, fp, "int");
//		ReadLabel(fp);
//		XGRead(&(thisgraph->time), sizeof(SCALAR), 1, fp, SCALAR_CHAR);
//		maxTime = max(maxTime, thisgraph->time);
//		minTime = min(minTime, thisgraph->time);
//		XGRead(&n, sizeof(int), 1, fp, "int");
//		XGRead(&m, sizeof(int), 1, fp, "int");
//		if (!x) {
//			x = (SCALAR*)malloc(m * sizeof(SCALAR));
//			y = (SCALAR*)malloc(n * sizeof(SCALAR));
//			z = (SCALAR**)malloc(m * sizeof(SCALAR*));
//			w = (SCALAR**)malloc(m * sizeof(SCALAR*));
//			for (int j = 0; j < m; j++) {
//				z[j] = (SCALAR*)calloc(n, sizeof(SCALAR));
//				w[j] = (SCALAR*)calloc(n, sizeof(SCALAR));
//			}
//
//		}
//		// read the x and y axes
//		XGRead(x, sizeof(SCALAR), m, fp, SCALAR_CHAR);
//		XGRead(y, sizeof(SCALAR), n, fp, SCALAR_CHAR);
//
//
//		thisgraph->z = (SCALAR**)malloc(m * sizeof(SCALAR*));
//		thisgraph->w = (SCALAR**)malloc(m * sizeof(SCALAR*));
//		for (int j = 0; j < m; j++) {
//			thisgraph->w[j] = (SCALAR*)malloc(n * sizeof(SCALAR));
//			thisgraph->z[j] = (SCALAR*)malloc(n * sizeof(SCALAR));
//			XGRead(thisgraph->w[j], sizeof(SCALAR), n, fp, SCALAR_CHAR);
//			XGRead(thisgraph->z[j], sizeof(SCALAR), n, fp, SCALAR_CHAR);
//		}
//		tmp.add(thisgraph);
//	}
//	fclose(fp);
//
//	graphdata = tmp;  //this will reverse tmp
//	current = new ListIter<VectorPlotData>(graphdata);
//	XGSetVec("vecvec", label->X_Label, label->Y_Label, label->Z_Label,
//		"open", xloc, yloc, label->X_Scale, label->Y_Scale,
//		label->X_Auto_Rescale, label->Y_Auto_Rescale,
//		label->X_Min, label->X_Max, label->Y_Min, label->Y_Max);
//	XGVector(x, y, w, z, &m, &n, 1);
//	ListIter<StructureData> walk2(structures);
//	for (walk2.restart(); !walk2.Done(); walk2++)
//		XGStructureArray(walk2.current()->n, (STRUCT_FILL)walk2.current()->fillFlag,
//			walk2.current()->lineColor,
//			walk2.current()->fillColor, walk2.current()->points);
//}

//void VectorPlot::updatePlot(double time) {
//
//	if (time == minTime) current->restart();
//
//	VectorPlotData* thisdata = current->current();
//	while (!current->Done() && current->current()->time <= time)
//	{
//		thisdata = current->current();
//		(*current)++;
//	}
//	if (thisdata)
//		for (int j = 0; j < m; j++) {
//			memcpy(z[j], thisdata->z[j], n * sizeof(SCALAR));
//			memcpy(w[j], thisdata->w[j], n * sizeof(SCALAR));
//		}
//}
//
//void XGMainLoop() {
//
//	//  cause the graphs to cycle
//	if (theTime > maxTime) theTime = minTime;
//
//	ListIter<plot> walk(*thePlots);
//	for (walk.restart(); !walk.Done(); walk++) {
//		walk.current()->updatePlot(theTime);
//	}
//	theTime += thetimestep;
//
//}

void Dump(const char*) {}
