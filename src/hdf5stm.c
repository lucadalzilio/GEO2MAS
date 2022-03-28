#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
//#undef _GNU_SOURCE

void create_output_hdf5( const char name[] )
{
	hid_t       file_id;   /* file identifier */
	herr_t      status;
	
	/* Create a new file using default properties. */
	file_id = H5Fcreate( name, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	
	//fprintf(fp_log,"Created: %s \n", name );
	
	/* Terminate access to the file. */
	status = H5Fclose(file_id);
}

void AddGroup_to_hdf5( const char filename[], const char group[] )
{
	hid_t       file_id, group_id;  /* identifiers */
	herr_t      status;
	char        *group_name;
	
	asprintf( &group_name, "/%s", group );
	
	/* Open exisiting file */
	file_id = H5Fopen( filename, H5F_ACC_RDWR, H5P_DEFAULT);
	
	/* Create group "ParticleGroup" in the root group using absolute name. */
#ifdef H5_1_8
	group_id = H5Gcreate(file_id, group_name, H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
#endif
#ifdef H5_1_6
	group_id = H5Gcreate(file_id, group_name, H5P_DEFAULT);
#endif
	
	/* Close group. */
	status = H5Gclose(group_id);
	
	/* Close the file. */
	status = H5Fclose(file_id);
	
	free( group_name );
}

void AddFieldToGroup_generic(
							 int compress,
							 const char filename[],
							 const char group[], const char field[],
							 char d_type,
							 int np, void* data, int dim )
{
	hid_t   file_id, particle_group_id, coord_dataset_id, coord_dataspace_id;  /* identifiers */
	herr_t  status;
	hsize_t length;
	int attr[2];
	hid_t attribute_id, attr_dataspace_id;
	hsize_t attr_dims;
	char *dataset_name;
	char *group_name;
	
	hid_t    plist=0,SET_CREATION_PLIST;
	hsize_t  cdims[2] = {0,0};
	double percentage_chunk;
	double chunk_size;
	int deflation_level;
	
	
	asprintf( &group_name, "/%s", group );
	asprintf( &dataset_name, "%s/%s", group, field );
	
	length = dim * np;
	
	/* Open exisiting file */
	file_id = H5Fopen( filename, H5F_ACC_RDWR, H5P_DEFAULT);
	
	/* Open group "ParticleGroup" */
#ifdef H5_1_8
	particle_group_id = H5Gopen(file_id, group_name,H5P_DEFAULT);
#endif
#ifdef H5_1_6
	particle_group_id = H5Gopen(file_id, group_name);
#endif
	
	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	
	percentage_chunk = 5.0;
	chunk_size = percentage_chunk*((double)np)/((double)100.0 );
	if( chunk_size < 1 ) {
		cdims[0] = 1;
	}
	else {
		cdims[0] = (int)chunk_size;
	}
	//      cdims[0] = 400;
	
	plist  = H5Pcreate(H5P_DATASET_CREATE);
	status = H5Pset_chunk(plist, 1, cdims);
	deflation_level = 4;
	status = H5Pset_deflate( plist, deflation_level);
	
	
	/* Create the data space for the dataset. */
	coord_dataspace_id = H5Screate_simple( 1, &length, NULL);
	
	if(compress==_TRUE_) {
		SET_CREATION_PLIST = plist;
		//fprintf(fp_log,"*** Compression info *** \n");
		//fprintf(fp_log,"  chunk_size = %f \n", chunk_size );
		//fprintf(fp_log,"  deflation level = %d \n", deflation_level );
	}
	else {
		SET_CREATION_PLIST = H5P_DEFAULT;
	}
	
	
	
	if( d_type == 'd' ) {
#ifdef H5_1_8
		/* Create a dataset within "group". */
		coord_dataset_id = H5Dcreate( file_id, dataset_name, H5T_NATIVE_DOUBLE, coord_dataspace_id, H5P_DEFAULT,SET_CREATION_PLIST,H5P_DEFAULT);
#endif
#ifdef H5_1_6
		coord_dataset_id = H5Dcreate( file_id, dataset_name, H5T_NATIVE_DOUBLE, coord_dataspace_id, SET_CREATION_PLIST );
#endif
		
		/* Write the particle dataset. */
		status = H5Dwrite(coord_dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data );
	}
	else if( d_type == 'i' ) {
#ifdef H5_1_8
		coord_dataset_id = H5Dcreate( file_id, dataset_name, H5T_STD_I32BE, coord_dataspace_id, H5P_DEFAULT,SET_CREATION_PLIST,H5P_DEFAULT);
#endif
#ifdef H5_1_6
		coord_dataset_id = H5Dcreate( file_id, dataset_name, H5T_STD_I32BE, coord_dataspace_id, SET_CREATION_PLIST);
#endif
		status = H5Dwrite(coord_dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data );
		
	}
	else if( d_type == 'f' ) {
#ifdef H5_1_8
		coord_dataset_id = H5Dcreate( file_id, dataset_name, H5T_NATIVE_FLOAT, coord_dataspace_id, H5P_DEFAULT,SET_CREATION_PLIST,H5P_DEFAULT);
#endif
#ifdef H5_1_6
		coord_dataset_id = H5Dcreate( file_id, dataset_name, H5T_NATIVE_FLOAT, coord_dataspace_id, SET_CREATION_PLIST);
#endif
		status = H5Dwrite(coord_dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data );
		
	}
	
	else {
		fprintf(fp_log,"ERROR: Only know how to write doubles (d), integers (i), and floats (f). \n");
		fflush(fp_log);
		exit(1);
	}
	/*Create attribute for particle field
	 "field_offset" tell us the stride to the next field associated with a given point. ie
	 data[ i*field_offset + d ], i is particle index, d is dof index (d<field_offset).*/
	attr_dims = 2;
	attr_dataspace_id = H5Screate_simple(1, &attr_dims, NULL);
	
	attr[0] = np;
	attr[1] = dim;
#ifdef H5_1_8
	attribute_id = H5Acreate(coord_dataset_id, "length,dim", H5T_STD_I32BE, attr_dataspace_id, H5P_DEFAULT,H5P_DEFAULT);
#endif
#ifdef H5_1_6
	attribute_id = H5Acreate(coord_dataset_id, "length,dim", H5T_STD_I32BE, attr_dataspace_id, H5P_DEFAULT);
#endif
	
	status = H5Awrite(attribute_id, H5T_NATIVE_INT, attr );
	status = H5Aclose(attribute_id);
	
	/* Close the data space for this particle dataset. */
	status = H5Sclose(attr_dataspace_id);
	
	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

	/* Close plist */
	status = H5Pclose(plist);
	
	/* Close the particle coord dataset. */
	status = H5Sclose(coord_dataspace_id);	
	
	/* Close the particle coord dataset. */
	status = H5Dclose(coord_dataset_id);
	
	/* Close group. */
	status = H5Gclose(particle_group_id);
	
	/* Close the file. */
	status = H5Fclose(file_id);
	
	free(group_name);
	free(dataset_name);
}

// Used for data storage during evolution to 'steady-state'; marker info is not saved
int create_hdf5( int compress, char* txtout, int n0 )
{
	double *mx, *my, arr[5], densgrid, chain_sz;
	//float *chain_x, *chain_y;
	int *ii, index, *compo_arr, opt, d;
	int  np=marknum, nn=nodenum, m1, ch_ind=0;
	char *name;
	
	/*compress = _FALSE_; // valgrind testing // */
	if(compress==_TRUE_) 
		{asprintf( &name, "%s%s",txtout,".gzip.h5");}
	else 
		{asprintf( &name, "%s%s" ,txtout, ".h5" );}
	create_output_hdf5( name );
	
 	AddGroup_to_hdf5( name, "ModelGroup" );
	AddGroup_to_hdf5( name, "NodeGroup" );
	// Interpolate marker properties on visualization grid
	AddGroup_to_hdf5( name, "VisMarkerGroup" );
	
	
	// Model Parameters
    	if (setup<10)
		{arr[0] = timesum;}
	else
		{arr[0] = timesum/3.15576e+7;;}
	arr[1] = xsize;
	arr[2] = ysize;
	arr[3] = xnumx;
	arr[4] = ynumy;
	
	// Composition on markers
	// Interpolates and retrieve composition on the dense marker grid (Taras)
	opt = 5;
	compo_arr = (int*)malloc( sizeof(int) * MAXMRK * 1 );
	memset(compo_arr, 0, MAXMRK*sizeof(int));
	// Compress rocktype = last 1
	index = interpoler(compo_arr, opt, n0, 1);
	
	mcomp = (int*)malloc( sizeof(int) * index * 1 );	
	for (m1=0;m1<index;m1++)
		{mcomp[m1] = compo_arr[m1];}
	
	// Put all the data in structured hdf5 file
	AddFieldToGroup_generic( compress, name, "ModelGroup" , "Model"   , 'd',     5,     arr, 1 );
	AddFieldToGroup_generic( compress, name, "ModelGroup" , "gx"      , 'd', xnumx,      gx, 1 );
	AddFieldToGroup_generic( compress, name, "ModelGroup" , "gy"      , 'd', ynumy,      gy, 1 );
	AddFieldToGroup_generic( compress, name, "VisMarkerGroup", "Mtype", 'i', index,   mcomp, 1 );
	AddFieldToGroup_generic( compress, name, "NodeGroup"  , "pr"      , 'd',    nn,      pr, 1 );
	AddFieldToGroup_generic( compress, name, "NodeGroup"  , "ro"      , 'd',    nn,      ro, 1 );
	AddFieldToGroup_generic( compress, name, "NodeGroup"  , "nu"      , 'd',    nn,      nu, 1 );
	AddFieldToGroup_generic( compress, name, "NodeGroup"  , "nd"      , 'd',    nn,      nd, 1 );
	AddFieldToGroup_generic( compress, name, "NodeGroup"  , "vx"      , 'd',    nn,      vx, 1 );
	AddFieldToGroup_generic( compress, name, "NodeGroup"  , "vy"      , 'd',    nn,      vy, 1 );
	AddFieldToGroup_generic( compress, name, "NodeGroup"  , "exx"     , 'd',    nn,     exx, 1 );
	AddFieldToGroup_generic( compress, name, "NodeGroup"  , "exy"     , 'd',    nn,     exy, 1 );
	AddFieldToGroup_generic( compress, name, "NodeGroup"  , "sxx"     , 'd',    nn,     sxx, 1 );
	AddFieldToGroup_generic( compress, name, "NodeGroup"  , "sxy"     , 'd',    nn,     sxy, 1 );
	AddFieldToGroup_generic( compress, name, "NodeGroup"  , "sbrit"   , 'd',    nn,  sbritn, 1 );
	AddFieldToGroup_generic( compress, name, "NodeGroup"  , "esp"     , 'd',    nn,     esp, 1 );
	AddFieldToGroup_generic( compress, name, "NodeGroup"  , "gg"      , 'd',    nn,      gg, 1 );
	AddFieldToGroup_generic( compress, name, "NodeGroup"  , "gd"      , 'd',    nn,      gd, 1 );
	if (setup>9)
		{AddFieldToGroup_generic( compress, name, "NodeGroup"  , "tk"      , 'd',    nn,      tk, 1 );}

	fprintf(fp_log, "File %s has been created\n", name );
	fflush(fp_log);
	free( compo_arr );
	free( name );
	free( mcomp );
	
	return 0;
}

// Store marker variables
int create_hdf5_markerprop( int compress, char* txtout, int n0 )
{
	int m1,m1f=0,m1ff=0,*tmp;
	char *name_me;
	double *tmpf;
	
	if(compress==_TRUE_) {
		asprintf( &name_me, "%s%s",txtout,".gzip.h5");
	}
	else {
		asprintf( &name_me, "%s%s" ,txtout, ".h5" );
	}
	create_output_hdf5( name_me );
	
	// Store individual marker properties
	AddGroup_to_hdf5( name_me, "MarkerGroup" );
	
	// --- Make array of which markers need to be saved; follow flag that sticks to marker ---
	// Only use those markers that were flagged as followable at dt = 1 of start_cond = 1
	// Markers that are later added, are added after marknum, so will automatically be added at the end of this new array

	// Resize memory if markers were added to array of interest (follow)
	if (nm > nm_old)
	{
		tmp = (int*)realloc(mcf, sizeof(int) * nm * 1 );
		if (tmp != NULL){
			mcf = tmp;
		} else {
			fprintf(stderr, "ERROR: Can't reallocate memory\n"); fflush(fp_log);
			// dynarray remains allocated 
			//mcf is still all there, so only but not overwrite mcf ,so do nothing else ?
			//free mcf q,     X store mcf = tmp q last can not because is dangling and unknown how much can contain
			//but can give other name ; can  never use that name anymore ??
			//can just alloc a bit bigger, and send nm with it in hdf5 ?
		}
		tmpf = (double*)realloc(markxf, sizeof(double) * nm * 1 );
		if (tmpf != NULL){
			markxf = tmpf;
		} else {
			fprintf(stderr, "ERROR: Can't reallocate memory\n"); fflush(fp_log);
			exit(EXIT_FAILURE);
		}
		tmpf = (double*)realloc(markyf, sizeof(double) * nm * 1 );
		if (tmpf != NULL){
			markyf = tmpf;
		} else {
			fprintf(stderr, "ERROR: Can't reallocate memory\n"); fflush(fp_log);
			exit(EXIT_FAILURE);
		}
		tmpf = (double*)realloc(markef, sizeof(double) * nm * 1 );
		if (tmpf != NULL){
			markef = tmpf;
		} else {
			// If have this problem think of something else, for suggestion see bottom or top this if. 
			fprintf(stderr, "ERROR: Can't reallocate memory\n"); fflush(fp_log);
			exit(EXIT_FAILURE);
		}
		tmpf = (double*)realloc(markexxf, sizeof(double) * nm * 1 );
		if (tmpf != NULL){
			markexxf = tmpf;
		} else {
			fprintf(stderr, "ERROR: Can't reallocate memory\n"); fflush(fp_log);
			exit(EXIT_FAILURE);
		}
		tmpf = (double*)realloc(markexyf, sizeof(double) * nm * 1 );
		if (tmpf != NULL){
			markexyf = tmpf;
		} else {
			fprintf(stderr, "ERROR: Can't reallocate memory\n"); fflush(fp_log);
			exit(EXIT_FAILURE);
		}
		tmpf = (double*)realloc(markxxf, sizeof(double) * nm * 1 );
		if (tmpf != NULL){
			markxxf = tmpf;
		} else {
			fprintf(stderr, "ERROR: Can't reallocate memory\n"); fflush(fp_log);
			exit(EXIT_FAILURE);
		}
		tmpf = (double*)realloc(markxyf, sizeof(double) * nm * 1 );
		if (tmpf != NULL){
			markxyf = tmpf;
		} else {
			fprintf(stderr, "ERROR: Can't reallocate memory\n"); fflush(fp_log);
			exit(EXIT_FAILURE);
		}
		tmpf = (double*)realloc(markvf, sizeof(double) * nm * 1 );
		if (tmpf != NULL){
			markvf = tmpf;
		} else {
			fprintf(stderr, "ERROR: Can't reallocate memory\n"); fflush(fp_log);
			exit(EXIT_FAILURE);
		}
		tmpf = (double*)realloc(markpf, sizeof(double) * nm * 1 );
		if (tmpf != NULL){
			markpf = tmpf;
		} else {
			fprintf(stderr, "ERROR: Can't reallocate memory\n"); fflush(fp_log);
			exit(EXIT_FAILURE);
		}
		nm_old=nm;
	}

	// Fill arrays for hdf5 storage	
	for (m1=0;m1<marknum;m1++)
	{
		if (follow[m1]==1)
		{
			// decide do not need these seperare t,x,y; since thrust thibault that vis comp works based on this save.c must work 
			mcf[m1f] = (int)(markt[m1]);
			markxf[m1f] = (markx[m1]);
			markyf[m1f] = (marky[m1]);
			markef[m1f] = (marke[m1]);
			markexxf[m1f] = (markexx[m1]);
			markexyf[m1f] = (markexy[m1]);
			markxxf[m1f] = (markxx[m1]);
			markxyf[m1f] = (markxy[m1]);
			markvf[m1f] = (markv[m1]);
			markpf[m1f] = (markp[m1]);
	
			m1f++;
		}
	}
	if (nm != m1f){ fprintf(fp_log,"WARNING: Number of values stored in hdf5-marker-arrays (%i) is not equal to number predicted (nm=%i) !!!\n",m1f,nm); fflush(fp_log); }	

	// Put all the data in structured hdf5 file
	AddFieldToGroup_generic( compress, name_me, "MarkerGroup", "markt"   , 'i',    nm,     mcf, 1 );
	AddFieldToGroup_generic( compress, name_me, "MarkerGroup", "markx"   , 'd',    nm,  markxf, 1 );
	AddFieldToGroup_generic( compress, name_me, "MarkerGroup", "marky"   , 'd',    nm,  markyf, 1 );
	AddFieldToGroup_generic( compress, name_me, "MarkerGroup", "marke"   , 'd',    nm,  markef, 1 );
	AddFieldToGroup_generic( compress, name_me, "MarkerGroup", "markexx" , 'd',    nm,markexxf, 1 );
	AddFieldToGroup_generic( compress, name_me, "MarkerGroup", "markexy" , 'd',    nm,markexyf, 1 );
	AddFieldToGroup_generic( compress, name_me, "MarkerGroup", "markxx"  , 'd',    nm, markxxf, 1 );
	AddFieldToGroup_generic( compress, name_me, "MarkerGroup", "markxy"  , 'd',    nm, markxyf, 1 );
	AddFieldToGroup_generic( compress, name_me, "MarkerGroup", "markp"   , 'd',    nm,  markpf, 1 );
	AddFieldToGroup_generic( compress, name_me, "MarkerGroup", "markv"   , 'd',    nm,  markvf, 1 );

	fprintf(fp_log, "File %s has been created (incl marker properties) \n", name_me );
	fflush(fp_log);
	free( name_me );
	
	return 0;
}

// Store fluid marker locations
int create_hdf5_fluid( int compress, char* txtout, int n0 )
{
	int m1,m1f=0,m1ff=0,*tmp;
	char *name_me;
	double *tmpf;
	
	if(compress==_TRUE_)
		{asprintf( &name_me, "%s%s",txtout,".gzip.h5");}
	else 
		{asprintf( &name_me, "%s%s" ,txtout, ".h5" );}
	create_output_hdf5( name_me );
	
	// Store individual marker properties
	AddGroup_to_hdf5( name_me, "FluidGroup" );
	
	// --- Make array of which markers need to be saved; follow flag that sticks to marker ---
	// Only use those markers that were flagged as followable at dt = 1 of start_cond = 1
	// Markers that are later added, are added after marknum, so will automatically be added at the end of this new array

	// Resize memory if markers were added to array of interest (follow)
	if (nmf > nmf_old)
	{
		tmp = (int*)realloc(mcff, sizeof(int) * nmf * 1 );
		if (tmp != NULL)
			{mcff = tmp;}
		else
			{fprintf(stderr, "ERROR: Can't reallocate memory\n"); fflush(fp_log);}
			// dynarray remains allocated 
			//mcf is still all there, so only but not overwrite mcf ,so do nothing else ?
			//free mcf q,     X store mcf = tmp q last can not because is dangling and unknown how much can contain
			//but can give other name ; can  never use that name anymore ??
			//can just alloc a bit bigger, and send nm with it in hdf5 ?
		
		tmpf = (double*)realloc(markxff, sizeof(double) * nmf * 1 );
		if (tmpf != NULL)
			{markxff = tmpf;} 
		else 
			{fprintf(stderr, "ERROR: Can't reallocate memory\n"); fflush(fp_log); exit(EXIT_FAILURE);}
		
		tmpf = (double*)realloc(markyff, sizeof(double) * nmf * 1 );
		if (tmpf != NULL)
			{markyff = tmpf;} 
		else 
			{fprintf(stderr, "ERROR: Can't reallocate memory\n"); fflush(fp_log); exit(EXIT_FAILURE);}
		
		tmpf = (double*)realloc(markwff, sizeof(double) * nmf * 1 );
		if (tmpf != NULL)
			{markwff = tmpf;} 
		else 
			{fprintf(stderr, "ERROR: Can't reallocate memory\n"); fflush(fp_log); exit(EXIT_FAILURE);}
		
		nmf_old=nmf;
	}

	// Fill arrays for hdf5 storage	
	for (m1=0;m1<marknum;m1++)
	{
		if (follow[m1]==2)
		{
			mcff[m1ff] = (int)(markt[m1]);
			markxff[m1ff] = (markx[m1]);
			markyff[m1ff] = (marky[m1]);
			markwff[m1ff] = (markw[m1]);
	
			m1ff++;
		}
	}
	if (nmf != m1ff){ fprintf(fp_log,"WARNING: Number of values stored in hdf5-marker-arrays (%i) is not equal to number predicted (nm=%i) !!!\n",m1ff,nmf); fflush(fp_log); }	

	// Put all the data in structured hdf5 file
	AddFieldToGroup_generic( compress, name_me, "FluidGroup", "fmarkt"   , 'i',    nmf,     mcff, 1 );
	AddFieldToGroup_generic( compress, name_me, "FluidGroup", "fmarkx"   , 'd',    nmf,  markxff, 1 );
	AddFieldToGroup_generic( compress, name_me, "FluidGroup", "fmarky"   , 'd',    nmf,  markyff, 1 );
	AddFieldToGroup_generic( compress, name_me, "FluidGroup", "fmarkw"   , 'd',    nmf,  markwff, 1 );
	
	fprintf(fp_log, "File %s has been created (fluid location) \n", name_me ); fflush(fp_log);
	free( name_me );
	
	return 0;
}

/* Print Results to data file ----------------------------------- */
int interpoler(int compo_arr[], int param, int n0, int compress_rt)
/* n0 - file number */
// param - which parameter to interpolate and output
{
	/* Counters */
	int n1,n2,mm2,mm3,ccur,cmin=0,cmax=0,yn,index=0,m10;
	char nn1;
	long int m1,m2,m3,m4,m5,m6,m7;
	long int mm1,erosmark=0;
	/* Buffers for XY */
	double x,y;
	char szint,szlong,szfloat,szdouble;
	double ival, ival1,xresx,yresy,xminx,xmaxx,yminy,ymaxy,ea,na,mpb,mtk,mro,mbb,mnu,mkt,mcp,mht,mvp,mvs,xcur,ycur,mwa,mdi0,maa,mgg,mxmelt,mhlatent;
	long int xresol,yresol,nresol,markstp;
	long int xresol0,xresolsum;
	double xincr,xminx0,yminy0;
	/* TD Database variables  */
	double n,e;
	/**/
	/**/

	//if (printmod) fprintf(fp_log,"Print data to %s...",fl1out);
	/**/
	/**/
	/**/
	/*
	 // Open file 
	 fl = fopen(fl1out,"wt");
	 */
	/* Read check resolution */
	/*
	 xincr=xresol=xresol0=(long int)(fl0stp[n0][2]);
	 yresol=(long int)(fl0stp[n0][5]);
	 nresol=xresol0*yresol;
	 */
	xincr=xsize/res_high+1;
	xresol0=xresol= (long int) xincr;
	yresol= (long int) (ysize/res_high+1);
	nresol=xresol0*yresol;
	/* Check coordinates, calc step *//*
									   if (fl0stp[n0][0]<0) fl0stp[n0][0]=0;
									   if (fl0stp[n0][1]>1.0) fl0stp[n0][1]=1.0;
									   if (fl0stp[n0][3]<0) fl0stp[n0][3]=0;
									   if (fl0stp[n0][4]>1.0) fl0stp[n0][4]=1.0;
									   xminx0=xminx=fl0stp[n0][0]*xsize;
									   xmaxx=fl0stp[n0][1]*xsize;
									   xresx=(xmaxx-xminx)/(double)(xresol-1);
									   yminy0=yminy=fl0stp[n0][3]*ysize;
									   ymaxy=fl0stp[n0][4]*ysize;
									   yresy=(ymaxy-yminy)/(double)(yresol-1);
									   yminy-=yresy; if(yminy<0) yminy=0;
									   ymaxy+=yresy; if(ymaxy>ysize) ymaxy=ysize;*/
	xminx0=xminx=0.0;
	xmaxx=xsize;
	/*xmaxx=fl0stp[n0][1]*xsize;*/
	xresx=(xmaxx-xminx)/(double)(xresol-1);
	/*yminy0=yminy=fl0stp[n0][3]*ysize;*/
	yminy0=yminy=0.0;
	/*ymaxy=fl0stp[n0][4]*ysize;*/
	ymaxy=ysize;
	yresy=(ymaxy-yminy)/(double)(yresol-1);
	yminy-=yresy; if(yminy<0) yminy=0;
	ymaxy+=yresy; if(ymaxy>ysize) ymaxy=ysize;
 
	if((nresol*2)>MAXMAT)
	{
		fprintf(fp_log,"\n Limited space in val0[] lin0[], spleating ...\n");
		fflush(fp_log);
		xincr=(long int)((double)(MAXMAT)/(double)(yresol*2));
	}
	/* X, Y Resolution save */
	/*
	 fprintf(fp_log,"%ld %ld %ld   %e %e %e %e %e %e\n",xresol,yresol,nresol,xminx,xmaxx,xresx,yminy,ymaxy,yresy); getchar();
	 */
	//fprintf(fl,"%e \n %ld %ld \n",timesum/3.15576e+7,xresol0,yresol);
	/* Step for markers definition */
	ival=0.1*(double)(marknum)/(double)(nresol);
	markstp=(long int)(ival);
	if (markstp<1) markstp=1;
	/**/
	/**/
	/**/
	/* Compute vorticity tensor */
	for (m1=1;m1<xnumx-1;m1++)
		for (m2=1;m2<ynumy-1;m2++)
		{
			/* Pos of in ol0[] */
			m3=m1*ynumy+m2;
			/**/
			/* Min,Max Vx definition */
			esp[m3]=spncalc(m1,m2);
			/*
			 fprintf(fp_log,"%d %ld %ld %ld %e %e %e %e ",param,m1,m2,m3,exy[m3],eps[0],eps[1],esp[m3]); getchar();
			 */
		}
	/* Compute vorticity tensor */
	/**/
	/**/
	/**/
	/* Save data in text format ---------------------------- */
	xresolsum=0;
	do
	{
		/* Set x resolution */
		xresol=xincr;
		if(xresolsum+xresol>xresol0) xresol=xresol0-xresolsum;
		nresol=xresol*yresol;
		//xminx=xminx0=fl0stp[n0][0]*xsize+xresx*((double)(xresolsum));
		xminx-=xresx; if(xminx<0) xminx=0;
		//xmaxx=fl0stp[n0][0]*xsize+xresx*((double)(xresolsum+xresol));
		if(xmaxx>xsize) xmaxx=xsize;
		/**/
		/* Clear visual arrays */
		for (m1=0;m1<nresol;m1++)
		{
			val0[m1]=val0[nresol+m1]=0;
			lin0[m1]=lin0[nresol+m1]=0;
		}
		/**/
		/**/
		/**/
		/* ROCK TYPE VISUALISATION ----------------------------*/
		/* CHEMICAL COMPONENT VISUALISATION ----------------------------*/
		/* Variations in markers types for chemical component visualisation */
		if(param==5) 
		{
			cmin=1000,cmax=0;
			/* Markers cycle */
			for (m3=0;m3<=marknum;m3+=markstp) 
			/* Check markers out of visualisation area */
				if ((markx[m3])>xminx && (markx[m3])<xmaxx && (marky[m3])>yminy && (marky[m3])<ymaxy && markt[m3]<50)
				{
					if (cmin>markt[m3]) cmin=markt[m3];
					if (cmax<markt[m3]) cmax=markt[m3];
				}
			/**/
			/* Marker type cycle for chemical component Visualisation */
			for (ccur=cmin;ccur<=cmax;ccur++) 
			{
				/* Markers cycle */
				for (m3=0;m3<=marknum;m3+=markstp) 
				/* Check markers out of visualisation area */
					if ((markx[m3])>xminx && (markx[m3])<xmaxx && (marky[m3])>yminy && (marky[m3])<ymaxy && markt[m3]==ccur)
					{
						/*
						 fprintf(fp_log,"%d %d %d   %ld %d   %e %e %e \n",cmin,cmax,ccur,m3,markt[m3],markx[m3],marky[m3],markk[m3]); getchar();
						 */
						/* Define relative position of the marker */
						ea=((markx[m3])-xminx0)/xresx+1.0;
						m1=(long int)(ea);
						if(m1<0) m1=0;
						if(m1>xresol) m1=xresol;
						ea-=(double)(m1);
						na=((marky[m3])-yminy0)/yresy+1.0;
						m2=(long int)(na);
						if(m2<0) m2=0;
						if(m2>yresol) m2=yresol;
						na-=(double)(m2);
						/* Add weights for four nodes surrounding the marker */
						m4=nresol+(m1-1)*yresol+(m2-1);
						if(m1>0 && m2>0) val0[m4]+=(float)((1.0-ea)*(1.0-na));
						if(m1>0 && m2<yresol) val0[m4+1]+=(float)((1.0-ea)*na);
						if(m1<xresol && m2>0) val0[m4+yresol]+=(float)(ea*(1.0-na));
						if(m1<xresol && m2<yresol) val0[m4+yresol+1]+=(float)(ea*na);
						/*
						 if(xresolsum){fprintf(fp_log,"%ld %ld %ld %ld  %e %e %e %e \n",m3,m1,m2,m4,markx[m3],marky[m3],ea,na); getchar();}
						 */
					}
				/**/
				/* Change types for nodes of visual arrays */
				for (m1=0;m1<nresol;m1++)
				{
					if(val0[nresol+m1]>val0[m1]) 
					{
						val0[m1]=val0[nresol+m1]; 
						lin0[m1]=ccur;
						/* Clear value */
						val0[nresol+m1]=0; 
					}
				}
			}
			/**/
			/* Marker type reload for chemical component Visualisation */
			yn=0;
			for (m1=0;m1<nresol;m1++)
			{
				if(!val0[m1]) 
				{
					lin0[m1]=-1;
					yn=1;
				}
			}
			/**/
			/* Marker type reinterpolate for empty nodes */
			if(gridmod && yn) 
			{
				/* Node cycle */
				for (m1=0;m1<xresol;m1++)
					for (m2=0;m2<yresol;m2++)
					{
						m3=m1*yresol+m2;
						/* Interpolate marker types from other nodes */
						if(lin0[m3]==-1) 
						{
							/* Search for surrounding non empty nodes */
							m4=1;
							yn=0;
							do
							{
								cmin=1000,cmax=0;
								for (m5=m1-m4;m5<=m1+m4;m5++)
									for (m6=m2-m4;m6<=m2+m4;m6++)
										if (m5>=0 && m5<xresol && m6>=0 && m6<yresol)
										{
											m7=m5*yresol+m6;
											if(lin0[m7]!=-1)
											{
												yn=1;
												if (cmin>lin0[m7]) cmin=lin0[m7];
												if (cmax<lin0[m7]) cmax=lin0[m7];
											}
										}
								m4++;
							}
							while(!yn && m4<=gridmod);
							/**/
							/* Recalc using non empty nodes */
							if(yn) 
							{
								/* Marker type cycle for chemical component Visualisation */
								for (ccur=cmin;ccur<=cmax;ccur++) 
								{
									for (m5=m1-m4+1;m5<m1+m4;m5++)
										for (m6=m2-m4+1;m6<m2+m4;m6++)
											if (m5>=0 && m5<xresol && m6>=0 && m6<yresol)
											{
												m7=m5*yresol+m6;
												/* Add weight */
												if(lin0[m7]==ccur)
												{
													ea=ABSV(((double)(m5-m1))/(double)(m4));
													na=ABSV(((double)(m6-m2))/(double)(m4));
													val0[nresol+m3]+=(float)((1.0-ea)*(1.0-na))*val0[m7];
													/*
													 val0[nresol+m3]+=(float)((1.0-ea)*(1.0-na));
													 {fprintf(fp_log,"%ld %ld %ld   %ld   %ld %ld %ld   %d \n",m1,m2,m3,m4,m5,m6,m7,ccur); getchar();}
													 */
												}
											}
									/* Set new marker type */
									if(val0[nresol+m3]>val0[m3]) 
									{
										val0[m3]=val0[nresol+m3]; 
										lin0[nresol+m3]=ccur;
										/* Clear value */
										val0[nresol+m3]=0; 
									}
								}
							}
						}
					}
				/**/
				/* Marker type reload for chemical component Visualisation */
				for (m1=0;m1<nresol;m1++)
				{
					if(lin0[m1]==-1 && val0[m1]) 
					{
						lin0[m1]=lin0[nresol+m1];
						lin0[nresol+m1]=0;
						yn=1;
					}
				}
			}
		}
		/* End CHEMICAL COMPONENT VISUALISATION ----------------------------*/
		/**/
		/**/
		/**/
		/* TRANSPORT PROPERTIES VISUALISATION +++++++++++++++++++++++++++++ */
		// If want to visualize marker properties relating to temperature, density, shear heating calculations ed, check large-scale code with routine from Thibault
		/* End TRANSPORT PROPERTIES VISUALISATION +++++++++++++++++++++++++++++ */
		/**/
		/**/
		/**/
		/* Visualisation parameters output */
		for (m1=0;m1<xresol;m1++)
			for (m2=0;m2<yresol;m2++)
			{
				/* Visual node number */
				m4=m1*yresol+m2;
				/* Visual node coordinates */
				xcur=xminx0+xresx*((double)(m1));
				ycur=yminy0+yresy*((double)(m2));
				/* Value for Visualisation */
				if(param==5 && compress_rt==1)
				{
					/*Compression of rock type */
					m5=1;
					// Check if similar rocktypes are aligned
					for (m3=m2+1;m3<yresol;m3++)
					{
						if(lin0[m4]==lin0[m4+m3-m2]) m5++; else break;
					}
					// If more than 3 next to each other; compress in 3 spaces give info on 1) if compressed (-2), 2) how long line, 3) rock type
					if(m5>3)
					{
						//fprintf(fl," %d %ld %ld",-2,m5,lin0[m4]);
						compo_arr[index] = -2;
						//fprintf(fp_log,"m5 is %d\n", m5);
						compo_arr[index+1] = m5;
						compo_arr[index+2] = (int)round(lin0[m4]);
						if(xresolsum < xresol-1){index+=3;}
						m2=m3-1;
					}
					// Not aligned, so not compressible, so just write 
					else
					{
						//fprintf(fl," %ld",lin0[m4]);
						compo_arr[index] = (int)round(lin0[m4]);
						if(xresolsum < xresol-1){index+=1;}
					}
				}
				else if(param==5 && compress_rt==0)
				{
					compo_arr[m4] = (int)round(lin0[m4]);
					//compo_arr[index] = (int)round(lin0[m4]);
					if(xresolsum < xresol-1){index+=1;}
				}
                		else
				{
					compo_arr[index] = (int)round(val0[m4]);
					if(xresolsum < xresol-1){index+=1;}
				}
			}
        /* Add x resolution */
        xresolsum+=xresol;
	}
	while(xresolsum<xresol0);
	return(index);
}


/* Spin tensor calc */ 
/* Spn=1/2(dVy/dX-dVx/dY) */
double spncalc(long int m1, long int m2)
/* m1,m2 - node X,Y number */
{
	/* Exy position */
	double xi,leftsxy=0,leftsxy1=0,leftsxy2=0;
	long int m1min,m1max,m2min,m2max,m3;
	int n1,nx,ny;
	/**/
	/**/
	/**/
	/* Staggered Nodes num */
	/*  [0]                [2]            */
	/*                                    */
	/*                     Vx2            */
	/*                                    */
	/*  [1]     Vy1        [3]       Vy3  */
	/*                   Exy3,Nu3         */
	/*                                    */
	/*                     Vx3            */
	/*                                    */
	/**/
	/**/
	/**/
	/* dVx/dY */
	xi=gy[m2];
	/* Calc, Check Fd limits */
	m2min=m2-1-stoksfd;
	if(m2min<0) m2min=0;
	m2max=m2+stoksfd;
	if(m2max>ynumy-2) m2max=ynumy-2;
	/**/
	/* Load distances to xn[] */
	for (m3=m2min;m3<=m2max;m3++)
	{
		xn[m3-m2min]=(gy[m3]+gy[m3+1])/2.0;
	}
	/**/
	/* Calc maximal position in xn[] */
	nx=(int)(m2max-m2min);
	/**/
	/* Calc Vx coefficients for EPSxy */
	fdweight(nx,1,xi);
	/* Reload coefficients to cn[] */
	for (m3=0;m3<=nx;m3++)
	{
		cn[m3][2]=cn[m3][1];
	}
	/**/
	/**/
	/**/
	/* dVy/dX */
	xi=gx[m1];
	/* Calc, Check Fd limits */
	m1min=m1-1-stoksfd;
	if(m1min<0) m1min=0;
	m1max=m1+stoksfd;
	if(m1max>xnumx-2) m1max=xnumx-2;
	/**/
	/* Load distances to xn[] */
	for (m3=m1min;m3<=m1max;m3++)
	{
		xn[m3-m1min]=(gx[m3]+gx[m3+1])/2.0;
	}
	/**/
	/* Calc maximal position in xn[] */
	ny=(int)(m1max-m1min);
	/**/
	/* Calc Vy coefficients for EPSxy */
	fdweight(ny,1,xi);
	/**/
	/* Return Spn val ----------------------------*/
	/* Exy=1/2(dVx/dY+dVy/dX)=0 */
	/* 1/2dVx/dY */
	/* Add Vx with koefficients */
	for (m3=m2min;m3<=m2max;m3++)
	{
		leftsxy1+=vx[m1*ynumy+m3]*cn[m3-m2min][2]/2.0;
	}
	/**/
	/* 1/2dVy/dX */
	/* Add Vy with koefficients */
	for (m3=m1min;m3<=m1max;m3++)
	{
		leftsxy2+=vy[m3*ynumy+m2]*cn[m3-m1min][1]/2.0;
	}
	/**/
	/* Save Exy */
	eps[0]=leftsxy=leftsxy1+leftsxy2;
	/* Save Esp (rotation rate) */
	eps[1]=leftsxy2-leftsxy1;

	return leftsxy2-leftsxy1;
}
/* Spin tensor calc */ 





