/* Load information from configuration file mode.t3c ============== */
int loadconf()
{
	/* Counter */
	int n1,n2,fln3=0;

	/* Open File file.t3c */
	fl = fopen("file.t3c","rt");
	ffscanf(); fln3=atoi(sa)-1;
	fclose(fl);

	/* Open File mode.t3c */
#if setup < 10
	fl = fopen("mode_lab.t3c","rt");
#elif setup > 9 && ls_t0 == 1
	fl = fopen("mode.t3c","rt");
#elif setup > 9 && ls_t0 == 0
	fl = fopen("mode_ls_stm.t3c","rt");
#endif
	
	/* Data File name */
	ffscanf();
	for (n1=0;n1<50;n1++) fl1in[n1]=sa[n1];
	ffscanf(); if(sa[0] == 'b') fl1itp=1;

	/* Load first Results File names */
	ffscanf();
	fl0num=0;
	while(sa[0]!='~')
	{
		/* Check file Counter */
		if(fl0num>=MAXFLN) {fprintf(fp_log,"Space out in fl0out[]"); fflush(fp_log); exit(0);}
		/**/
		/* Save results file name */
		for (n1=0;n1<50;n1++) fl0out[fl0num][n1]=sa[n1];
		/**/
		/* Load TYPE_cyc0max_maxxystep_maxtkstep_maxtmstep_maxdsstep */
		ffscanf(); if(sa[0] == 'b') {fl0otp[fl0num]=1;} 
		else if(sa[0] == 'h') {fl0otp[fl0num]=2;}
		ffscanf();fl0cyc[fl0num]=atoi(sa);
		ffscanf();fl0stp[fl0num][0]=atof(sa);
		ffscanf();fl0stp[fl0num][1]=atof(sa);
		ffscanf();fl0stp[fl0num][2]=atof(sa)*3.15576e+7;
		ffscanf();fl0stp[fl0num][3]=atof(sa)*3.15576e+7;
		/**/
		/* Incr File Counters */
		fl0num++;
		/**/
		/* Load Next Results File names */
		ffscanf();
	}
	/**/
	/* Data File name change after number */
	if(fln3>=0 && fln3<fl0num)
	{
		for (n1=0;n1<50;n1++) fl1in[n1]=fl0out[fln3][n1];
		fl1itp=fl0otp[fln3];
	}
	else
	{
		fln3=-1;
	}
	/**/
	/* Service */
	ffscanf();printmod=atoi(sa);
	ffscanf();timedir=atof(sa);
	ffscanf();movemod=atoi(sa);
	ffscanf();tempmod=atoi(sa);
	ffscanf();markmod=atoi(sa);
	ffscanf();ratemod=atoi(sa);
	ffscanf();gridmod=atoi(sa);
	ffscanf();intermod=atoi(sa);
	ffscanf();intermod1=atoi(sa);
	ffscanf();outgrid=atoi(sa);
	ffscanf();densimod=atoi(sa);
        ffscanf();timebond=atof(sa)*3.15576e+7;
	/**/
	/* Erosion/Sedimentation */
	ffscanf();erosmod=atoi(sa);
	ffscanf();eroslev=atof(sa);
	ffscanf();eroscon=atof(sa);
	ffscanf();eroskoe=atof(sa);
	ffscanf();sedilev=atof(sa);
	ffscanf();sedicon=atof(sa);
	ffscanf();sedikoe=atof(sa);
	ffscanf();sedimcyc=atoi(sa);
	ffscanf();waterlev=atof(sa);
	ffscanf();slopemax=atof(sa);
	/**/
	/* V */
	ffscanf();DIVVMIN=atof(sa);
	ffscanf();STOKSMIN=atof(sa);
	ffscanf();stoksmod=atoi(sa);
	ffscanf();viscmod=atof(sa);
	ffscanf();stoksfd=atoi(sa);
	ffscanf();nubeg=atof(sa);
	ffscanf();nuend=atof(sa);
	ffscanf();nucontr=atof(sa);
	ffscanf();hidry=atof(sa);
	ffscanf();hidrl=atof(sa);
	ffscanf();strmin=atof(sa);
	ffscanf();strmax=atof(sa);
	ffscanf();stredif=atof(sa);
	/**/
	/**/
	/**/
	/* T */
	ffscanf();HEATMIN=atof(sa);
	ffscanf();heatmod=atoi(sa);
	ffscanf();heatfd=atoi(sa);
	ffscanf();heatdif=atof(sa);
	ffscanf();frictyn=atoi(sa);
	ffscanf();adiabyn=atoi(sa);
	/**/
	/**/
	/**/
	/* Water */
	ffscanf();tkpor=atof(sa);
	ffscanf();zmpor=atof(sa);
	ffscanf();vyfluid=atof(sa);
	ffscanf();vymelt=atof(sa);
	ffscanf();dmwamin=atof(sa);
	ffscanf();tdeep=atof(sa);
	ffscanf();zdeep=atof(sa);
	ffscanf();dxwater=atof(sa);
	ffscanf();dywater=atof(sa);
	ffscanf();deserp=atof(sa);
	ffscanf();dyserp=atof(sa);
	/**/
	/**/
	/**/
	/* Seismo-Thermo-Mechanical parameters */
	ffscanf();pemod=atoi(sa);
	ffscanf();savefluid=atof(sa);
	ffscanf();meltmod=atof(sa);
	ffscanf();debugmod=atof(sa);
	/* eclogitization - temp range */	
 	ffscanf();teclmin=atof(sa);
	ffscanf();teclmax=atof(sa);
        ffscanf();collision=atof(sa);
	/*
	fprintf(fp_log,"%e %e %e %e %e",tkpor,zmpor,vyfluid,vymelt,dmwamin);getchar();
	*/
	/**/
	fclose(fl);
	/* End Load information from configuration file mode.t3c */


	/* stop.yn file creation */
	fl = fopen("stop.yn","wt");
	fprintf(fl,"n \n");
	fclose(fl);
	/**/
	//return fln3;


        if (densimod>=2)
        {

	/* Load thermodynamic databases */
#if setup>9
	if (printmod){fprintf(fp_log,"Loading thermodynamic databases... \n"); fflush(fp_log);}
		
	/* Dry peridotite */
	/* RO - density */
	fl = fopen("./Database/pdry_rho","rt");
	ffscanf();
	ffscanf(); tknum=atoi(sa);
	ffscanf(); pbnum=atoi(sa);
	ffscanf(); tkmin=atof(sa);
	ffscanf(); pbmin=atof(sa);
	ffscanf(); tkstp=atof(sa);
	ffscanf(); pbstp=atof(sa);
	ffscanf();
	ffscanf();
	//	fprintf(fp_log,"%d %d %e %e %e %e",tknum,pbnum,tkmin,pbmin,tkstp,pbstp);getchar();
	
	for (n1=0;n1<pbnum;n1++)
		for (n2=0;n2<tknum;n2++)
	{
		ffscanf(); td[n2][n1][0][0]=atof(sa)/1000.0;
	}
	fclose(fl);

	/* H  - enthalpy */
	fl = fopen("./Database/pdry_hh","rt");
	for (n1=0;n1<9;n1++) ffscanf();
	for (n1=0;n1<pbnum;n1++)
		for (n2=0;n2<tknum;n2++)
	{
		ffscanf(); td[n2][n1][0][1]=atof(sa)/1000.0/4.1837;
	}
	fclose(fl);

	/* Vp - seismic velosity, km/s */
	fl = fopen("./Database/pdry_vp","rt");
	for (n1=0;n1<9;n1++) ffscanf();
	for (n1=0;n1<pbnum;n1++)
		for (n2=0;n2<tknum;n2++)
	{
		ffscanf(); td[n2][n1][0][2]=MAXV(0,atof(sa));
	}
	fclose(fl);

	/* Vs - seismic velosity, km/s */
	fl = fopen("./Database/pdry_vs","rt");
	for (n1=0;n1<9;n1++) ffscanf();
	for (n1=0;n1<pbnum;n1++)
		for (n2=0;n2<tknum;n2++)
	{
		ffscanf(); td[n2][n1][0][3]=MAXV(0,atof(sa));
	}
	fclose(fl);
	
	fprintf(fp_log,"pdry_vs OK \n");

	/* Wet peridotite */
	/* RO - density */
	fl = fopen("./Database/pwet_rho","rt");
	for (n1=0;n1<9;n1++) ffscanf();
	for (n1=0;n1<pbnum;n1++)
		for (n2=0;n2<tknum;n2++)
	{
		ffscanf(); td[n2][n1][1][0]=atof(sa)/1000.0;
	}
	fclose(fl);
	/**/
	/* H  - enthalpy */
	fl = fopen("./Database/pwet_hh","rt");
	for (n1=0;n1<9;n1++) ffscanf();
	for (n1=0;n1<pbnum;n1++)
		for (n2=0;n2<tknum;n2++)
	{
		ffscanf(); td[n2][n1][1][1]=atof(sa)/1000.0/4.1837;
	}
	fclose(fl);
	/**/
	/* Vp - seismic velosity, km/s */
	fl = fopen("./Database/pwet_vp","rt");
	for (n1=0;n1<9;n1++) ffscanf();
	for (n1=0;n1<pbnum;n1++)
		for (n2=0;n2<tknum;n2++)
	{
		ffscanf(); td[n2][n1][1][2]=MAXV(0,atof(sa));
	}
	fclose(fl);

	/* Vs - seismic velosity, km/s */
	fl = fopen("./Database/pwet_vs","rt");
	for (n1=0;n1<9;n1++) ffscanf();
	for (n1=0;n1<pbnum;n1++)
		for (n2=0;n2<tknum;n2++)
	{
		ffscanf(); td[n2][n1][1][3]=MAXV(0,atof(sa));
	}
	fclose(fl);

	/* Wa - water contents, wt% */
	fl = fopen("./Database/pwet_h2o","rt");
	for (n1=0;n1<9;n1++) ffscanf();
	for (n1=0;n1<pbnum;n1++)
		for (n2=0;n2<tknum;n2++)
	{
		ffscanf(); td[n2][n1][1][4]=MAXV(0,atof(sa));
	}
	fclose(fl);
	fprintf(fp_log,"pwet_h2o OK \n");

	/* Molten peridotite */
	/* RO - density */
	fl = fopen("./Database/pwetmelt_rho","rt");
	for (n1=0;n1<9;n1++) ffscanf();
	for (n1=0;n1<pbnum;n1++)
		for (n2=0;n2<tknum;n2++)
	{
		ffscanf(); td[n2][n1][2][0]=atof(sa)/1000.0;
	}
	fclose(fl);

	/* H  - enthalpy */
	fl = fopen("./Database/pwetmelt_hh","rt");
	for (n1=0;n1<9;n1++) ffscanf();
	for (n1=0;n1<pbnum;n1++)
		for (n2=0;n2<tknum;n2++)
	{
		ffscanf(); td[n2][n1][2][1]=atof(sa)/1000.0/4.1837;
	}
	fclose(fl);
	/**/
	/* Vp - seismic velosity, km/s */
	fl = fopen("./Database/pwetmelt_vp","rt");
	for (n1=0;n1<9;n1++) ffscanf();
	for (n1=0;n1<pbnum;n1++)
		for (n2=0;n2<tknum;n2++)
	{
		ffscanf(); td[n2][n1][2][2]=MAXV(0,atof(sa));
	}
	fclose(fl);

	/* Vs - seismic velosity, km/s */
	fl = fopen("./Database/pwetmelt_vs","rt");
	for (n1=0;n1<9;n1++) ffscanf();
	for (n1=0;n1<pbnum;n1++)
		for (n2=0;n2<tknum;n2++)
	{
		ffscanf(); td[n2][n1][2][3]=MAXV(0,atof(sa));
	}
	fclose(fl);

	/* Wa - water contents, wt% */
	fl = fopen("./Database/pwetmelt_h2o","rt");
	for (n1=0;n1<9;n1++) ffscanf();
	for (n1=0;n1<pbnum;n1++)
		for (n2=0;n2<tknum;n2++)
	{
		ffscanf(); td[n2][n1][2][4]=MAXV(0,atof(sa));
	}
	fclose(fl);
	fprintf(fp_log,"pwetmelt_h2o OK \n");


	/* Wet Gabbro */
	/* RO - density */
	fl = fopen("./Database/gabwet_rho","rt");
	for (n1=0;n1<9;n1++) ffscanf();
	for (n1=0;n1<pbnum;n1++)
		for (n2=0;n2<tknum;n2++)
	{
		ffscanf(); td[n2][n1][3][0]=atof(sa)/1000.0;
	}
	fclose(fl);

	/* H  - enthalpy */
	fl = fopen("./Database/gabwet_hh","rt");
	for (n1=0;n1<9;n1++) ffscanf();
	for (n1=0;n1<pbnum;n1++)
		for (n2=0;n2<tknum;n2++)
	{
		ffscanf(); td[n2][n1][3][1]=atof(sa)/1000.0/4.1837;
		/*
		if (n2>0) 
		{
		ival=1000.0*4.1837*(td[n2][n1][3][1]-td[n2-1][n1][3][1])/tkstp;
		if (ival<=1e+2 || ival>=5e+4) {fprintf(fp_log,"GWET %d %d %e <%s>",n1,n2,ival,sa);getchar();}
		}
		*/
	}
	fclose(fl);

	/* Vp - seismic velosity, km/s */
	fl = fopen("./Database/gabwet_vp","rt");
	for (n1=0;n1<9;n1++) ffscanf();
	for (n1=0;n1<pbnum;n1++)
		for (n2=0;n2<tknum;n2++)
	{
		ffscanf(); td[n2][n1][3][2]=MAXV(0,atof(sa));
	}
	fclose(fl);

	/* Vs - seismic velosity, km/s */
	fl = fopen("./Database/gabwet_vs","rt");
	for (n1=0;n1<9;n1++) ffscanf();
	for (n1=0;n1<pbnum;n1++)
		for (n2=0;n2<tknum;n2++)
	{
		ffscanf(); td[n2][n1][3][3]=MAXV(0,atof(sa));
	}
	fclose(fl);

	/* Wa - water contents, wt% */
	fl = fopen("./Database/gabwet_h2o","rt");
	for (n1=0;n1<9;n1++) ffscanf();
	for (n1=0;n1<pbnum;n1++)
		for (n2=0;n2<tknum;n2++)
	{
		ffscanf(); td[n2][n1][3][4]=MAXV(0,atof(sa));
	}
	fprintf(fp_log,"gabwet_h2o OK \n");
	fclose(fl);


	/* Molten Gabbro */
	/* RO - density */
	fl = fopen("./Database/gabwetmelt_rho","rt");
	for (n1=0;n1<9;n1++) ffscanf();
	for (n1=0;n1<pbnum;n1++)
		for (n2=0;n2<tknum;n2++)
	{
		ffscanf(); td[n2][n1][4][0]=atof(sa)/1000.0;
	}
	fclose(fl);

	/* H  - enthalpy */
	fl = fopen("./Database/gabwetmelt_hh","rt");
	for (n1=0;n1<9;n1++) ffscanf();
	for (n1=0;n1<pbnum;n1++)
		for (n2=0;n2<tknum;n2++)
	{
		ffscanf(); td[n2][n1][4][1]=atof(sa)/1000.0/4.1837;
	}
	fclose(fl);

	/* Vp - seismic velosity, km/s */
	fl = fopen("./Database/gabwetmelt_vp","rt");
	for (n1=0;n1<9;n1++) ffscanf();
	for (n1=0;n1<pbnum;n1++)
		for (n2=0;n2<tknum;n2++)
	{
		ffscanf(); td[n2][n1][4][2]=MAXV(0,atof(sa));
	}
	fclose(fl);

	/* Vs - seismic velosity, km/s */
	fl = fopen("./Database/gabwetmelt_vs","rt");
	for (n1=0;n1<9;n1++) ffscanf();
	for (n1=0;n1<pbnum;n1++)
		for (n2=0;n2<tknum;n2++)
	{
		ffscanf(); td[n2][n1][4][3]=MAXV(0,atof(sa));
	}
	fclose(fl);

	/* Wa - water contents, wt% */
	fl = fopen("./Database/gabwetmelt_h2o","rt");
	for (n1=0;n1<9;n1++) ffscanf();
	for (n1=0;n1<pbnum;n1++)
		for (n2=0;n2<tknum;n2++)
	{
		ffscanf(); td[n2][n1][4][4]=MAXV(0,atof(sa));
	}
	fclose(fl);
	fprintf(fp_log,"gabwetmelt_h2o OK \n");


	/* Wet sediments */
	/* RO - density */
	fl = fopen("./Database/swet_rho","rt");
	for (n1=0;n1<9;n1++) ffscanf();
	for (n1=0;n1<pbnum;n1++)
		for (n2=0;n2<tknum;n2++)
	{
		ffscanf(); td[n2][n1][5][0]=atof(sa)/1000.0;
	}
	fclose(fl);

	/* H  - enthalpy */
	fl = fopen("./Database/swet_hh","rt");
	for (n1=0;n1<9;n1++) ffscanf();
	for (n1=0;n1<pbnum;n1++)
		for (n2=0;n2<tknum;n2++)
	{
		ffscanf(); td[n2][n1][5][1]=atof(sa)/1000.0/4.1837;
	}
	fclose(fl);

	/* Vp - seismic velosity, km/s */
	fl = fopen("./Database/swet_vp","rt");
	for (n1=0;n1<9;n1++) ffscanf();
	for (n1=0;n1<pbnum;n1++)
		for (n2=0;n2<tknum;n2++)
	{
		ffscanf(); td[n2][n1][5][2]=MAXV(0,atof(sa));
	}
	fclose(fl);

	/* Vs - seismic velosity, km/s */
	fl = fopen("./Database/swet_vs","rt");
	for (n1=0;n1<9;n1++) ffscanf();
	for (n1=0;n1<pbnum;n1++)
		for (n2=0;n2<tknum;n2++)
	{
		ffscanf(); td[n2][n1][5][3]=MAXV(0,atof(sa));
	}
	fclose(fl);

	/* Wa - water contents, wt% */
	fl = fopen("./Database/swet_h2o","rt");
	for (n1=0;n1<9;n1++) ffscanf();
	for (n1=0;n1<pbnum;n1++)
		for (n2=0;n2<tknum;n2++)
	{
		ffscanf(); td[n2][n1][5][4]=MAXV(0,atof(sa));
	}
	fclose(fl);
	fprintf(fp_log,"swet_h2o OK \n");


	/* Molten sediments */
	/* RO - density */
	fl = fopen("./Database/swetmelt_rho","rt");
	for (n1=0;n1<9;n1++) ffscanf();
	for (n1=0;n1<pbnum;n1++)
		for (n2=0;n2<tknum;n2++)
	{
		ffscanf(); td[n2][n1][6][0]=atof(sa)/1000.0;
	}
	fclose(fl);
	/**/
	/* H  - enthalpy */
	fl = fopen("./Database/swetmelt_hh","rt");
	for (n1=0;n1<9;n1++) ffscanf();
	for (n1=0;n1<pbnum;n1++)
		for (n2=0;n2<tknum;n2++)
	{
		ffscanf(); td[n2][n1][6][1]=atof(sa)/1000.0/4.1837;
	}
	fclose(fl);

	/* Vp - seismic velosity, km/s */
	fl = fopen("./Database/swetmelt_vp","rt");
	for (n1=0;n1<9;n1++) ffscanf();
	for (n1=0;n1<pbnum;n1++)
		for (n2=0;n2<tknum;n2++)
	{
		ffscanf(); td[n2][n1][6][2]=MAXV(0,atof(sa));
	}
	fclose(fl);

	/* Vs - seismic velosity, km/s */
	fl = fopen("./Database/swetmelt_vs","rt");
	for (n1=0;n1<9;n1++) ffscanf();
	for (n1=0;n1<pbnum;n1++)
		for (n2=0;n2<tknum;n2++)
	{
		ffscanf(); td[n2][n1][6][3]=MAXV(0,atof(sa));
	}
	fclose(fl);

	/* Wa - water contents, wt% */
	fl = fopen("./Database/swetmelt_h2o","rt");
	for (n1=0;n1<9;n1++) ffscanf();
	for (n1=0;n1<pbnum;n1++)
		for (n2=0;n2<tknum;n2++)
	{
		ffscanf(); td[n2][n1][6][4]=MAXV(0,atof(sa));
	}
	fclose(fl);
	fprintf(fp_log,"swetmelt_h2o OK \n");
	
	
	/* Wet Basalt */
	/* RO - density */
	fl = fopen("./Database/bwet_rho","rt");
	for (n1=0;n1<9;n1++) ffscanf();
	for (n1=0;n1<pbnum;n1++)
		for (n2=0;n2<tknum;n2++)
	{
		ffscanf(); td[n2][n1][7][0]=atof(sa)/1000.0;
	}
	fclose(fl);

	/* H  - enthalpy */
	fl = fopen("./Database/bwet_hh","rt");
	for (n1=0;n1<9;n1++) ffscanf();
	for (n1=0;n1<pbnum;n1++)
		for (n2=0;n2<tknum;n2++)
	{
		ffscanf(); td[n2][n1][7][1]=atof(sa)/1000.0/4.1837;
	}
	fclose(fl);

	/* Vp - seismic velosity, km/s */
	fl = fopen("./Database/bwet_vp","rt");
	for (n1=0;n1<9;n1++) ffscanf();
	for (n1=0;n1<pbnum;n1++)
		for (n2=0;n2<tknum;n2++)
	{
		ffscanf(); td[n2][n1][7][2]=MAXV(0,atof(sa));
	}
	fclose(fl);

	/* Vs - seismic velosity, km/s */
	fl = fopen("./Database/bwet_vs","rt");
	for (n1=0;n1<9;n1++) ffscanf();
	for (n1=0;n1<pbnum;n1++)
		for (n2=0;n2<tknum;n2++)
	{
		ffscanf(); td[n2][n1][7][3]=MAXV(0,atof(sa));
	}
	fclose(fl);

	/* Wa - water contents, wt% */
	fl = fopen("./Database/bwet_h2o","rt");
	for (n1=0;n1<9;n1++) ffscanf();
	for (n1=0;n1<pbnum;n1++)
		for (n2=0;n2<tknum;n2++)
	{
		ffscanf(); td[n2][n1][7][4]=MAXV(0,atof(sa));
	}
	fclose(fl);
	fprintf(fp_log,"bwet_h2o OK \n");


	/* Molten Basalt */
	/* RO - density */
	fl = fopen("./Database/bwetmelt_rho","rt");
	for (n1=0;n1<9;n1++) ffscanf();
	for (n1=0;n1<pbnum;n1++)
		for (n2=0;n2<tknum;n2++)
	{
		ffscanf(); td[n2][n1][8][0]=atof(sa)/1000.0;
	}
	fclose(fl);

	/* H  - enthalpy */
	fl = fopen("./Database/bwetmelt_hh","rt");
	for (n1=0;n1<9;n1++) ffscanf();
	for (n1=0;n1<pbnum;n1++)
		for (n2=0;n2<tknum;n2++)
	{
		ffscanf(); td[n2][n1][8][1]=atof(sa)/1000.0/4.1837;
	}
	fclose(fl);

	/* Vp - seismic velosity, km/s */
	fl = fopen("./Database/bwetmelt_vp","rt");
	for (n1=0;n1<9;n1++) ffscanf();
	for (n1=0;n1<pbnum;n1++)
		for (n2=0;n2<tknum;n2++)
	{
		ffscanf(); td[n2][n1][8][2]=MAXV(0,atof(sa));
	}
	fclose(fl);

	/* Vs - seismic velosity, km/s */
	fl = fopen("./Database/bwetmelt_vs","rt");
	for (n1=0;n1<9;n1++) ffscanf();
	for (n1=0;n1<pbnum;n1++)
		for (n2=0;n2<tknum;n2++)
	{
		ffscanf(); td[n2][n1][8][3]=MAXV(0,atof(sa));
	}
	fclose(fl);

	/* Wa - water contents, wt% */
	fl = fopen("./Database/bwetmelt_h2o","rt");
	for (n1=0;n1<9;n1++) ffscanf();
	for (n1=0;n1<pbnum;n1++)
		for (n2=0;n2<tknum;n2++)
	{
		ffscanf(); td[n2][n1][8][4]=MAXV(0,atof(sa));
	}
	fclose(fl);
	fprintf(fp_log,"bwetmelt_h2o OK \n");


	/* Dry Upper crust */
	/* RO - density */
	fl = fopen("./Database/ucdry_rho","rt");
	for (n1=0;n1<9;n1++) ffscanf();
	for (n1=0;n1<pbnum;n1++)
		for (n2=0;n2<tknum;n2++)
	{
		ffscanf(); td[n2][n1][11][0]=atof(sa)/1000.0;
	}
	fclose(fl);

	/* H  - enthalpy */
	fl = fopen("./Database/ucdry_hh","rt");
	for (n1=0;n1<9;n1++) ffscanf();
	for (n1=0;n1<pbnum;n1++)
		for (n2=0;n2<tknum;n2++)
	{
		ffscanf(); td[n2][n1][11][1]=atof(sa)/1000.0/4.1837;
	}
	fclose(fl);

	/* Vp - seismic velosity, km/s */
	fl = fopen("./Database/ucdry_vp","rt");
	for (n1=0;n1<9;n1++) ffscanf();
	for (n1=0;n1<pbnum;n1++)
		for (n2=0;n2<tknum;n2++)
	{
		ffscanf(); td[n2][n1][8][2]=MAXV(0,atof(sa));
	}
	fclose(fl);

	/* Vs - seismic velosity, km/s */
	fl = fopen("./Database/ucdry_vs","rt");
	for (n1=0;n1<9;n1++) ffscanf();
	for (n1=0;n1<pbnum;n1++)
		for (n2=0;n2<tknum;n2++)
	{
		ffscanf(); td[n2][n1][11][3]=MAXV(0,atof(sa));
	}
	fclose(fl);
	fprintf(fp_log,"ucdry_vs OK \n");

	/* Wet Upper crust */
	/* RO - density */
	fl = fopen("./Database/ucwet_rho","rt");
	for (n1=0;n1<9;n1++) ffscanf();
	for (n1=0;n1<pbnum;n1++)
		for (n2=0;n2<tknum;n2++)
	{
		ffscanf(); td[n2][n1][12][0]=atof(sa)/1000.0;
	}
	fclose(fl);

	/* H  - enthalpy */
	fl = fopen("./Database/ucwet_hh","rt");
	for (n1=0;n1<9;n1++) ffscanf();
	for (n1=0;n1<pbnum;n1++)
		for (n2=0;n2<tknum;n2++)
	{
		ffscanf(); td[n2][n1][12][1]=atof(sa)/1000.0/4.1837;
	}
	fclose(fl);

	/* Vp - seismic velosity, km/s */
	fl = fopen("./Database/ucwet_vp","rt");
	for (n1=0;n1<9;n1++) ffscanf();
	for (n1=0;n1<pbnum;n1++)
		for (n2=0;n2<tknum;n2++)
	{
		ffscanf(); td[n2][n1][12][2]=MAXV(0,atof(sa));
	}
	fclose(fl);

	/* Vs - seismic velosity, km/s */
	fl = fopen("./Database/ucwet_vs","rt");
	for (n1=0;n1<9;n1++) ffscanf();
	for (n1=0;n1<pbnum;n1++)
		for (n2=0;n2<tknum;n2++)
	{
		ffscanf(); td[n2][n1][12][3]=MAXV(0,atof(sa));
	}
	fclose(fl);

	/* Wa - water contents, wt% */
	fl = fopen("./Database/ucwet_h2o","rt");
	for (n1=0;n1<9;n1++) ffscanf();
	for (n1=0;n1<pbnum;n1++)
		for (n2=0;n2<tknum;n2++)
	{
		ffscanf(); td[n2][n1][12][4]=MAXV(0,atof(sa));
	}
	fclose(fl);
	fprintf(fp_log,"ucwet_h2o OK \n");

	/* Dry Lower crust */
	/* RO - density */
	fl = fopen("./Database/lcdry_rho","rt");
	for (n1=0;n1<9;n1++) ffscanf();
	for (n1=0;n1<pbnum;n1++)
		for (n2=0;n2<tknum;n2++)
	{
		ffscanf(); td[n2][n1][13][0]=atof(sa)/1000.0;
	}
	fclose(fl);

	/* H  - enthalpy */
	fl = fopen("./Database/lcdry_hh","rt");
	for (n1=0;n1<9;n1++) ffscanf();
	for (n1=0;n1<pbnum;n1++)
		for (n2=0;n2<tknum;n2++)
	{
		ffscanf(); td[n2][n1][13][1]=atof(sa)/1000.0/4.1837;
	}
	fclose(fl);

	/* Vp - seismic velosity, km/s */
	fl = fopen("./Database/lcdry_vp","rt");
	for (n1=0;n1<9;n1++) ffscanf();
	for (n1=0;n1<pbnum;n1++)
		for (n2=0;n2<tknum;n2++)
	{
		ffscanf(); td[n2][n1][13][2]=MAXV(0,atof(sa));
	}
	fclose(fl);

	/* Vs - seismic velosity, km/s */
	fl = fopen("./Database/lcdry_vs","rt");
	for (n1=0;n1<9;n1++) ffscanf();
	for (n1=0;n1<pbnum;n1++)
		for (n2=0;n2<tknum;n2++)
	{
		ffscanf(); td[n2][n1][13][3]=MAXV(0,atof(sa));
	}
	fclose(fl);
	fprintf(fp_log,"lcdry_vs OK \n");

	/* Wet Lower crust */
	/* RO - density */
	fl = fopen("./Database/lcwet_rho","rt");
	for (n1=0;n1<9;n1++) ffscanf();
	for (n1=0;n1<pbnum;n1++)
		for (n2=0;n2<tknum;n2++)
	{
		ffscanf(); td[n2][n1][14][0]=atof(sa)/1000.0;
	}
	fclose(fl);

	/* H  - enthalpy */
	fl = fopen("./Database/lcwet_hh","rt");
	for (n1=0;n1<9;n1++) ffscanf();
	for (n1=0;n1<pbnum;n1++)
		for (n2=0;n2<tknum;n2++)
	{
		ffscanf(); td[n2][n1][14][1]=atof(sa)/1000.0/4.1837;
	}
	fclose(fl);

	/* Vp - seismic velosity, km/s */
	fl = fopen("./Database/lcwet_vp","rt");
	for (n1=0;n1<9;n1++) ffscanf();
	for (n1=0;n1<pbnum;n1++)
		for (n2=0;n2<tknum;n2++)
	{
		ffscanf(); td[n2][n1][14][2]=MAXV(0,atof(sa));
	}
	fclose(fl);

	/* Vs - seismic velosity, km/s */
	fl = fopen("./Database/lcwet_vs","rt");
	for (n1=0;n1<9;n1++) ffscanf();
	for (n1=0;n1<pbnum;n1++)
		for (n2=0;n2<tknum;n2++)
	{
		ffscanf(); td[n2][n1][14][3]=MAXV(0,atof(sa));
	}
	fclose(fl);

	/* Wa - water contents, wt% */
	fl = fopen("./Database/lcwet_h2o","rt");
	for (n1=0;n1<9;n1++) ffscanf();
	for (n1=0;n1<pbnum;n1++)
		for (n2=0;n2<tknum;n2++)
	{
		ffscanf(); td[n2][n1][14][4]=MAXV(0,atof(sa));
	}
	fclose(fl);
	fprintf(fp_log,"lcwet_h2o OK \n");



	/* Load global thermodynamic database */
	/* Dry peridotite Xmg=0.895 */
	/* RO - density */
	fl = fopen("./Database/m895_ro","rt");
	ffscanf();
	ffscanf(); tknum1=atoi(sa);
	ffscanf(); pbnum1=atoi(sa);
	ffscanf(); tkmin1=atof(sa);
	ffscanf(); pbmin1=atof(sa);
	ffscanf(); tkstp1=atof(sa);
	ffscanf(); pbstp1=atof(sa);
	ffscanf();
	ffscanf();
	for (n1=0;n1<pbnum1;n1++)
		for (n2=0;n2<tknum1;n2++)
	{
		ffscanf(); td[n2][n1][9][0]=atof(sa)/1000.0;
	}
	fclose(fl);

	/* H  - enthalpy */
	fl = fopen("./Database/m895_hh","rt");
	for (n1=0;n1<9;n1++) ffscanf();
	for (n1=0;n1<pbnum1;n1++)
		for (n2=0;n2<tknum1;n2++)
	{
		ffscanf(); td[n2][n1][9][1]=atof(sa)/1000.0/4.1837;
	}
	fclose(fl);

	/* Vp - seismic velosity, km/s */
	fl = fopen("./Database/m895_vp","rt");
	for (n1=0;n1<9;n1++) ffscanf();
	for (n1=0;n1<pbnum1;n1++)
		for (n2=0;n2<tknum1;n2++)
	{
		ffscanf(); td[n2][n1][9][2]=MAXV(0,atof(sa));
	}
	fclose(fl);

	/* Vs - seismic velosity, km/s */
	fl = fopen("./Database/m895_vs","rt");
	for (n1=0;n1<9;n1++) ffscanf();
	for (n1=0;n1<pbnum1;n1++)
		for (n2=0;n2<tknum1;n2++)
	{
		ffscanf(); td[n2][n1][9][3]=MAXV(0,atof(sa));
	}
	fclose(fl);
	fprintf(fp_log,"m895_vs OK \n");

	/* Wet peridotite */
	/* RO - density */
	fl = fopen("./Database/morn_ro","rt");
	for (n1=0;n1<9;n1++) ffscanf();
	for (n1=0;n1<pbnum1;n1++)
		for (n2=0;n2<tknum1;n2++)
	{
		ffscanf(); td[n2][n1][10][0]=atof(sa)/1000.0;
	}
	fclose(fl);

	/* H  - enthalpy */
	fl = fopen("./Database/morn_hh","rt");
	for (n1=0;n1<9;n1++) ffscanf();
	for (n1=0;n1<pbnum1;n1++)
		for (n2=0;n2<tknum1;n2++)
	{
		ffscanf(); td[n2][n1][10][1]=atof(sa)/1000.0/4.1837;
	}
	fclose(fl);

	/* Vp - seismic velosity, km/s */
	fl = fopen("./Database/morn_vp","rt");
	for (n1=0;n1<9;n1++) ffscanf();
	for (n1=0;n1<pbnum1;n1++)
		for (n2=0;n2<tknum1;n2++)
	{
		ffscanf(); td[n2][n1][10][2]=MAXV(0,atof(sa));
	}
	fclose(fl);

	/* Vs - seismic velosity, km/s */
	fl = fopen("./Database/morn_vs","rt");
	for (n1=0;n1<9;n1++) ffscanf();
	for (n1=0;n1<pbnum1;n1++)
		for (n2=0;n2<tknum1;n2++)
	{
		ffscanf(); td[n2][n1][10][3]=MAXV(0,atof(sa));
	}
	fclose(fl);
	fprintf(fp_log,"morn_vs OK \n");
	fflush(fp_log);
	
	if (printmod){fprintf(fp_log,"Finished loading thermodynamic databases... \n"); fflush(fp_log);}
	}
#endif
	return fln3;
}


/* Load information from prn file ============== */
void loader()
	/* bondv[] - bondary value */
	/* bondm[] - bondary mode 0=Not, -1=Value, 1,2...=LinNum+1 */
	/* m1,m2 - node X,Y number */
{
	/* Counter */
	int n1,ir,im;
	char nn1;
	long int m1,m2,m3;
	long int mm1;
	char szint,szlong,szfloat,szdouble,szcur;
	double ival0,ival1;
	float ivalf;
	FILE *flb;

	/* Load Past Results from data file-------------------------------- */
	if (printmod) fprintf(fp_log,"Load past results from %s ...\n",fl1in); fflush(fp_log);

	/* Load in Binary Format ---------------------------- */
	flb = fopen(fl1in,"rb");
	/**/
	/* Sizes of var definition */
	szint=sizeof(n1);
	szlong=sizeof(m1);
	szfloat=sizeof(ivalf);
	szdouble=sizeof(ival1);
	
	/* Check sizes of variables */
	fread(&szcur,1,1,flb);
	if (szcur!=szint) {fprintf(fp_log,"Current INT size <%d> is different from given in file <%d> \n",szint,szcur); fflush(fp_log); exit(0);}
	fread(&szcur,1,1,flb);
	if (szcur!=szlong) {fprintf(fp_log,"Current LONG INT size <%d> is different from given in file <%d> \n",szlong,szcur); fflush(fp_log); exit(0);}
	fread(&szcur,1,1,flb);
	if (szcur!=szfloat) {fprintf(fp_log,"Current FLOAT size <%d> is different from given in file <%d> \n",szfloat,szcur); fflush(fp_log); exit(0);}
	fread(&szcur,1,1,flb);
	if (szcur!=szdouble) {fprintf(fp_log,"Current DOUBLE size <%d> is different from given in file <%d> \n",szdouble,szcur); fflush(fp_log); exit(0);}

	/* Grid Parameters */
	fread(&xnumx,szlong,1,flb);
	fread(&ynumy,szlong,1,flb);
	fread(&mnumx,szlong,1,flb);
	fread(&mnumy,szlong,1,flb);
	fread(&marknum,szlong,1,flb);
	fread(&xsize,szdouble,1,flb);
	fread(&ysize,szdouble,1,flb);
	fread(&pinit,szdouble,1,flb);
	fread(pkf,szdouble,4,flb);
	fread(&GXKOEF,szdouble,1,flb);
	fread(&GYKOEF,szdouble,1,flb);
	fread(&rocknum,szint,1,flb);
	fread(&bondnum,szlong,1,flb);
	fread(&n1,szint,1,flb);
	fread(&timesum,szdouble,1,flb);timesum*=3.15576e+7;

	/* Calc,Check Grid parameters */
	gridcheck();

	/* Rock Types information */
	fread(markim,szint,rocknum,flb);
	fread(markn0,szdouble,rocknum,flb);
	fread(markn1,szdouble,rocknum,flb);
	fread(marks0,szdouble,rocknum,flb);
	fread(marks1,szdouble,rocknum,flb);
	fread(marknu,szdouble,rocknum,flb);
	fread(markdh,szdouble,rocknum,flb);
	fread(markdv,szdouble,rocknum,flb);
	fread(markss,szdouble,rocknum,flb);
	fread(markmm,szdouble,rocknum,flb);
	fread(markgg,szdouble,rocknum,flb);
	fread(markll,szdouble,rocknum,flb);
	fread(marka0,szdouble,rocknum,flb);
	fread(marka1,szdouble,rocknum,flb);
	fread(markb0,szdouble,rocknum,flb);
	fread(markb1,szdouble,rocknum,flb);
	fread(marke0,szdouble,rocknum,flb);
	fread(marke1,szdouble,rocknum,flb);
	fread(markf0,szdouble,rocknum,flb);
	fread(markf1,szdouble,rocknum,flb);
	fread(markro,szdouble,rocknum,flb);
	fread(markbb,szdouble,rocknum,flb);
	fread(markaa,szdouble,rocknum,flb);
	fread(markcp,szdouble,rocknum,flb);
	fread(markkt,szdouble,rocknum,flb);
	fread(markkf,szdouble,rocknum,flb);
	fread(markkp,szdouble,rocknum,flb);
	fread(markht,szdouble,rocknum,flb);

	/* Nodes information */
	/* Vx,Vy,bondm[],ro[],nu[],ep[],et[],tk[],cp[],kt[],ht[] */
	for (m1=0;m1<nodenum;m1++)
	{
		fread(&ival0,szdouble,1,flb);pr[m1]=(ival0);
		fread(&ival0,szdouble,1,flb);vx[m1]=(ival0);
		fread(&ival0,szdouble,1,flb);vy[m1]=(ival0);
		fread(&m2,szlong,1,flb);bondm[m1*3+0]=m2;
		fread(&m2,szlong,1,flb);bondm[m1*3+1]=m2;
		fread(&m2,szlong,1,flb);bondm[m1*3+2]=m2;
		fread(&ival0,szdouble,1,flb);exx[m1]=(ival0);
		fread(&ival0,szdouble,1,flb);dro[m1]=(ival0);
		fread(&ival0,szdouble,1,flb);exy[m1]=(ival0);
		fread(&ival0,szdouble,1,flb);sxx[m1]=(ival0);
		fread(&ival0,szdouble,1,flb);drp[m1]=(ival0);
		fread(&ival0,szdouble,1,flb);sxy[m1]=(ival0);
		fread(&ival0,szdouble,1,flb);ro[m1]=(ival0);
		fread(&ival0,szdouble,1,flb);nu[m1]=(ival0);
		fread(&ival0,szdouble,1,flb);nd[m1]=(ival0);
		fread(&ival0,szdouble,1,flb);sxxe[m1]=(ival0);
		fread(&ival0,szdouble,1,flb);sppe[m1]=(ival0);
		fread(&ival0,szdouble,1,flb);sxye[m1]=(ival0);
		fread(&ival0,szdouble,1,flb);exxe[m1]=(ival0);
		fread(&ival0,szdouble,1,flb);exye[m1]=(ival0);
		fread(&ival0,szdouble,1,flb);esp[m1]=(ival0);
		fread(&ival0,szdouble,1,flb);mu[m1]=(ival0);
		fread(&ival0,szdouble,1,flb);gg[m1]=(ival0);
		fread(&ival0,szdouble,1,flb);gd[m1]=(ival0);
		fread(&ival0,szdouble,1,flb);mvx[m1]=(ival0);
		fread(&ival0,szdouble,1,flb);mvy[m1]=(ival0);
		fread(&ival0,szdouble,1,flb);mrx[m1]=(ival0);
		fread(&ival0,szdouble,1,flb);mry[m1]=(ival0);
		fread(&ival0,szdouble,1,flb);ep[m1]=(ival0);
		fread(&ival0,szdouble,1,flb);et[m1]=(ival0);
		fread(&ival0,szdouble,1,flb);tk[m1]=(ival0);
		fread(&m2,szlong,1,flb);bondm[nodenum3+m1]=m2;
		fread(&ival0,szdouble,1,flb);cp[m1]=(ival0);
		fread(&ival0,szdouble,1,flb);kt[m1]=(ival0);
		fread(&ival0,szdouble,1,flb);ht[m1]=(ival0);
	}

	/* Gridlines positions */
	// Load high resolution area properties for usage latter on
#if (setup<9)
	res_high = 100.0e3;
	fread(&ival0,szdouble,1,flb);gx[0]=(ival0);
	for (m1=1;m1<xnumx;m1++)
	{
		fread(&ival0,szdouble,1,flb);gx[m1]=(ival0);
		// Get minumum grid size / maximum resolution
		if (gx[m1]-gx[m1-1]<res_high)
			{ res_high = gx[m1]-gx[m1-1]; }
	}
	for (m2=0;m2<ynumy;m2++)
		{fread(&ival0,szdouble,1,flb);gy[m2]=(ival0);}
#else
	// m10_hr - m11_hr
	// m20_hr
	//   |
	// m21_hr               
	res_high = 100.0e3;
	for (m1=0;m1<xnumx;m1++)
	{
		fread(&ival0,szdouble,1,flb);gx[m1]=(ival0);
		// Get High Resolution area nodal limits
		if (m1>1)
		{
			if (gx[m1]-gx[m1-1]==gx[m1-1]-gx[m1-2] && gx[m1]-gx[m1-1]<res_high)
			{
				m10_hr = m1-2;
				res_high = gx[m1]-gx[m1-1];
			}
			if (gx[m1]-gx[m1-1]>gx[m1-1]-gx[m1-2] && gx[m1-1]-gx[m1-2]==res_high)
			{
				m11_hr = m1-1;
			}
		}
	}
	//fprintf(fp_log,"Final X HR lim:  %e %i %i ", res_high, m10_hr, m11_hr); getchar();
	res_high = 100.0e3;
	for (m2=0;m2<ynumy;m2++)
	{
		fread(&ival0,szdouble,1,flb);gy[m2]=(ival0);
		// Get High Resolution area nodal limits
		if (m2>1)
		{
			if (gy[m2]-gy[m2-1]==gy[m2-1]-gy[m2-2] && gy[m2]-gy[m2-1]<res_high)
			{
				m20_hr = m2-2;
				res_high = gy[m2]-gy[m2-1];
			}
			if (gy[m2]-gy[m2-1]>gy[m2-1]-gy[m2-2] && gy[m2-1]-gy[m2-2]==res_high)
			{
				m21_hr = m2-1;
			}
		}
	}
	//fprintf(fp_log,"Final Y HR lim:  %e %i %i ", res_high, m20_hr, m21_hr); getchar();
#endif

	/* Boundary Conditions Equations */
	for (m1=1;m1<bondnum;m1++)
	{
		fread(&ival0,szdouble,1,flb);bondv[m1][0]=(ival0);
		fread(&ival0,szdouble,1,flb);bondv[m1][1]=(ival0);
		fread(&ival0,szdouble,1,flb);bondv[m1][2]=(ival0);
		fread(&ival0,szdouble,1,flb);bondv[m1][3]=(ival0);
		fread(&m2,szlong,1,flb);bondn[m1][0]=m2;
		fread(&m2,szlong,1,flb);bondn[m1][1]=m2;
		fread(&m2,szlong,1,flb);bondn[m1][2]=m2;
	}

	/* Markers X,Y,types */
	for (mm1=0;mm1<marknum;mm1++)
	{
		fread(&ival0,szdouble,1,flb);markx[mm1]=ival0;
		fread(&ival0,szdouble,1,flb);marky[mm1]=ival0;
		fread(&ival0,szdouble,1,flb);markk[mm1]=ival0;
		fread(&ival0,szdouble,1,flb);marke[mm1]=ival0;
		fread(&ival0,szdouble,1,flb);markxx[mm1]=ival0;
		fread(&ival0,szdouble,1,flb);markv[mm1]=ival0;
		fread(&ival0,szdouble,1,flb);markxy[mm1]=ival0;
		fread(&ival0,szdouble,1,flb);markp[mm1]=ival0;
		fread(&ival0,szdouble,1,flb);markexx[mm1]=ival0;
		fread(&ival0,szdouble,1,flb);markexy[mm1]=ival0;
		fread(&ival0,szdouble,1,flb);markd[mm1]=ival0;
		fread(&ival0,szdouble,1,flb);markw[mm1]=ival0;
		fread(&ival0,szdouble,1,flb);markvx[mm1]=ival0;
		fread(&ival0,szdouble,1,flb);markvy[mm1]=ival0;
		fread(&nn1,1,1,flb);markt[mm1]=nn1;
	}
	
	fclose(flb);
	/* End: Load in Binary Format ---------------------------- */
	
	// Load additional data that can potentially be changed on the fly	
	// & Select markers that will be tracked each timestep
#if setup<10
	flyinput_lab();
#elif setup>9
	flyinput_dynw();
#endif	

	if (printmod) fprintf(fp_log,"Loading from prn-file OK!\n"); fflush(fp_log);
	if (debugmod) fprintf(fp_log,"Loading ok ? gx(1000)=%e markx(1000)=%e markp(1000)=%e \n",gx[1000],markx[1000],markp[1000]); fflush(fp_log);
}


/* Save information to prn file ------------------------------- */
void saver(int f0, int n0)
	/* n0 - circle number */
{
	/* Counters */
	int n1,mm2;
	char nn1;
	long int m1,m2,m3;
	long int mm1;
	/* Buffers */
	char szint,szlong,szfloat,szdouble;
	double ival0,ival1;
	float ivalf;

	if (printmod) fprintf(fp_log,"Print %d circle results to %s...",n0+1,fl1out); fflush(fp_log);	

	/* Save data in binary format ---------------------------- */
	fl = fopen(fl1out,"wb");

	/* Sizes of var definition */
	szint=sizeof(n1);
	szlong=sizeof(m1);
	szfloat=sizeof(ivalf);
	szdouble=sizeof(ival0);
	fwrite(&szint,1,1,fl);
	fwrite(&szlong,1,1,fl);
	fwrite(&szfloat,1,1,fl);
	fwrite(&szdouble,1,1,fl);

	/* Grid Parameters */
	fwrite(&xnumx,szlong,1,fl);
	fwrite(&ynumy,szlong,1,fl);
	fwrite(&mnumx,szlong,1,fl);
	fwrite(&mnumy,szlong,1,fl);
	fwrite(&marknum,szlong,1,fl);
	fwrite(&xsize,szdouble,1,fl);
	fwrite(&ysize,szdouble,1,fl);
	fwrite(&pinit,szdouble,1,fl);
	fwrite(pkf,szdouble,4,fl);
	fwrite(&GXKOEF,szdouble,1,fl);
	fwrite(&GYKOEF,szdouble,1,fl);
	fwrite(&rocknum,szint,1,fl);
	fwrite(&bondnum,szlong,1,fl);
	fwrite(&n0,szint,1,fl);
	ival1=timesum/3.15576e+7;fwrite(&ival1,szdouble,1,fl);

	/* Rock Types information */
	fwrite(markim,szint,rocknum,fl);
	fwrite(markn0,szdouble,rocknum,fl);
	fwrite(markn1,szdouble,rocknum,fl);
	fwrite(marks0,szdouble,rocknum,fl);
	fwrite(marks1,szdouble,rocknum,fl);
	fwrite(marknu,szdouble,rocknum,fl);
	fwrite(markdh,szdouble,rocknum,fl);
	fwrite(markdv,szdouble,rocknum,fl);
	fwrite(markss,szdouble,rocknum,fl);
	fwrite(markmm,szdouble,rocknum,fl);
	fwrite(markgg,szdouble,rocknum,fl);
	fwrite(markll,szdouble,rocknum,fl);
	fwrite(marka0,szdouble,rocknum,fl);
	fwrite(marka1,szdouble,rocknum,fl);
	fwrite(markb0,szdouble,rocknum,fl);
	fwrite(markb1,szdouble,rocknum,fl);
	fwrite(marke0,szdouble,rocknum,fl);
	fwrite(marke1,szdouble,rocknum,fl);
	fwrite(markf0,szdouble,rocknum,fl);
	fwrite(markf1,szdouble,rocknum,fl);
	fwrite(markro,szdouble,rocknum,fl);
	fwrite(markbb,szdouble,rocknum,fl);
	fwrite(markaa,szdouble,rocknum,fl);
	fwrite(markcp,szdouble,rocknum,fl);
	fwrite(markkt,szdouble,rocknum,fl);
	fwrite(markkf,szdouble,rocknum,fl);
	fwrite(markkp,szdouble,rocknum,fl);
	fwrite(markht,szdouble,rocknum,fl);

	/* Nodes information */
	/* Vx,Vy,bondm[],ro[],nu[],ep[],et[],tk[],cp[],kt[],ht[] */
	for (m1=0;m1<nodenum;m1++)
	{
		ival0=(pr[m1]);fwrite(&ival0,szdouble,1,fl);
		ival0=(vx[m1]);fwrite(&ival0,szdouble,1,fl);
		ival0=(vy[m1]);fwrite(&ival0,szdouble,1,fl);
		m2=bondm[m1*3+0];fwrite(&m2,szlong,1,fl);
		m2=bondm[m1*3+1];fwrite(&m2,szlong,1,fl);
		m2=bondm[m1*3+2];fwrite(&m2,szlong,1,fl);
		ival0=(exx[m1]);fwrite(&ival0,szdouble,1,fl);
		ival0=(dro[m1]);fwrite(&ival0,szdouble,1,fl);
		ival0=(exy[m1]);fwrite(&ival0,szdouble,1,fl);
		ival0=(sxx[m1]);fwrite(&ival0,szdouble,1,fl);
		ival0=(drp[m1]);fwrite(&ival0,szdouble,1,fl);
		ival0=(sxy[m1]);fwrite(&ival0,szdouble,1,fl);
		ival0=(ro[m1]);fwrite(&ival0,szdouble,1,fl);
		ival0=(nu[m1]);fwrite(&ival0,szdouble,1,fl);
		ival0=(nd[m1]);fwrite(&ival0,szdouble,1,fl);
		ival0=(sxxe[m1]);fwrite(&ival0,szdouble,1,fl);
		ival0=(sppe[m1]);fwrite(&ival0,szdouble,1,fl);
		ival0=(sxye[m1]);fwrite(&ival0,szdouble,1,fl);
		ival0=(exxe[m1]);fwrite(&ival0,szdouble,1,fl);
		ival0=(exye[m1]);fwrite(&ival0,szdouble,1,fl);
		ival0=(esp[m1]);fwrite(&ival0,szdouble,1,fl);
		ival0=(mu[m1]);fwrite(&ival0,szdouble,1,fl);
		ival0=(gg[m1]);fwrite(&ival0,szdouble,1,fl);
		ival0=(gd[m1]);fwrite(&ival0,szdouble,1,fl);
		ival0=(mvx[m1]);fwrite(&ival0,szdouble,1,fl);
		ival0=(mvy[m1]);fwrite(&ival0,szdouble,1,fl);
		ival0=(mrx[m1]);fwrite(&ival0,szdouble,1,fl);
		ival0=(mry[m1]);fwrite(&ival0,szdouble,1,fl);
		ival0=(ep[m1]);fwrite(&ival0,szdouble,1,fl);
		ival0=(et[m1]);fwrite(&ival0,szdouble,1,fl);
		ival0=(tk[m1]);fwrite(&ival0,szdouble,1,fl);
		m2=bondm[nodenum3+m1];fwrite(&m2,szlong,1,fl);
		ival0=(cp[m1]);fwrite(&ival0,szdouble,1,fl);
		ival0=(kt[m1]);fwrite(&ival0,szdouble,1,fl);
		ival0=(ht[m1]);fwrite(&ival0,szdouble,1,fl);
	}

	/* Gridlines positions */
	for (m1=0;m1<xnumx;m1++)
	{
		ival0=(gx[m1]);fwrite(&ival0,szdouble,1,fl);
	}
	for (m2=0;m2<ynumy;m2++)
	{
		ival0=(gy[m2]);fwrite(&ival0,szdouble,1,fl);
	}

	/* Boundary Conditions Equations */
	for (m1=1;m1<bondnum;m1++)
	{
		ival0=(bondv[m1][0]);fwrite(&ival0,szdouble,1,fl);
		ival0=(bondv[m1][1]);fwrite(&ival0,szdouble,1,fl);
		ival0=(bondv[m1][2]);fwrite(&ival0,szdouble,1,fl);
		ival0=(bondv[m1][3]);fwrite(&ival0,szdouble,1,fl);
		m2=bondn[m1][0];fwrite(&m2,szlong,1,fl);
		m2=bondn[m1][1];fwrite(&m2,szlong,1,fl);
		m2=bondn[m1][2];fwrite(&m2,szlong,1,fl);
	}

	/* Markers X,Y,types */
	for (mm1=0;mm1<marknum;mm1++)
	{
		ival0=markx[mm1];fwrite(&ival0,szdouble,1,fl);
		ival0=marky[mm1];fwrite(&ival0,szdouble,1,fl);
		ival0=markk[mm1];fwrite(&ival0,szdouble,1,fl);
		ival0=marke[mm1];fwrite(&ival0,szdouble,1,fl);
		ival0=markxx[mm1];fwrite(&ival0,szdouble,1,fl);
		ival0=markv[mm1];fwrite(&ival0,szdouble,1,fl);
		ival0=markxy[mm1];fwrite(&ival0,szdouble,1,fl);
		ival0=markp[mm1];fwrite(&ival0,szdouble,1,fl);
		ival0=markexx[mm1];fwrite(&ival0,szdouble,1,fl);
		ival0=markexy[mm1];fwrite(&ival0,szdouble,1,fl);
		ival0=markd[mm1];fwrite(&ival0,szdouble,1,fl);
		ival0=markw[mm1];fwrite(&ival0,szdouble,1,fl);
		ival0=markvx[mm1];fwrite(&ival0,szdouble,1,fl);
		ival0=markvy[mm1];fwrite(&ival0,szdouble,1,fl);
		nn1=markt[mm1];fwrite(&nn1,1,1,fl);
	}

	fclose(fl);
	if (printmod) fprintf(fp_log,"Saved OK!\n"); fflush(fp_log);



	/* file.t3c file creation */
	fl = fopen("file.t3c","wt");
	fprintf(fl,"%d \n",f0);
	fclose(fl);



	/* stop.yn file information read */
	fl = fopen("stop.yn","rt");
	/* Read String */
	ffscanf();
	/* Stop Y/N */
	if (sa[0]=='y' || sa[0]=='Y')
	{
		fclose(fl);
		fprintf(fp_log,"PROGRAM TERMINATED FROM stop.yn \n");
		fflush(fp_log);	
		exit(0);
	}

	
	
	/* Change printmod */
	if (sa[0]>='0' && sa[0]<='9')
	{
		printmod=atoi(sa);
	}
	fclose(fl);
}


/* LOAD WITHOUT EMPTY LINES from fl =================================== */
void ffscanf()
{
	/* Counter */
	int n1;
	/**/
	/* Read cycle */
	do
	{
		/* Input string */
		n1=fscanf(fl,"%s",sa);
		
		/* Check end of file */
		if (n1<1)
		{
			fprintf(fp_log,"\n Unexpected end of file: while reading using ffscanf\n");
			fflush(fp_log);
			fclose(fl);
			exit(0);
		}
		
		/* Delete last symbol <32 */
		for(n1=strlen(sa)-1;n1>=0;n1--)
			if (*(sa+n1)<=32)
				*(sa+n1)=0;
		else
			break;
	}
	while (strlen(sa)<1 || sa[0]=='/');
}
/* End LOAD WITHOUT EMPTY LINES from fl =================================== */


/* LOAD WITHOUT EMPTY LINES from fl1 =================================== */
void ffscanf1()
{
	/* Counter */
	int n1;
	/**/
	/* Read cycle */
	do
	{
		/* Input string */
		n1=fscanf(fl1,"%s",sa);
		
		/* Check end of file */
		if (n1<1)
		{
			fprintf(fp_log,"\n Unexpected end of file: while reading using ffscanf1 (double)\n");
			fflush(fp_log);	
			fclose(fl1);
			exit(0);
		}
		
		/* Delete last symbol <32 */
		for(n1=strlen(sa)-1;n1>=0;n1--)
			if (*(sa+n1)<=32)
				*(sa+n1)=0;
		else
			break;
	}
	while (strlen(sa)<1 || sa[0]=='/');
}
/* End LOAD WITHOUT EMPTY LINES from fl1 =================================== */


/* Calc,Check Parameters of Grid */
void gridcheck()
{
	/* Gridlines NUM */
	if(xnumx>MAXXLN) {fprintf(fp_log,"Space out in gx[] %ld",xnumx); fflush(fp_log); exit(0);}
	if(ynumy>MAXYLN) {fprintf(fp_log,"Space out in gy[] %ld",ynumy); fflush(fp_log); exit(0);}
	/**/
	/* Nodes Num */
	nodenum=xnumx*ynumy;
	if(nodenum>MAXNOD) {fprintf(fp_log,"Space out in vx[],vy[] %ld",nodenum); fflush(fp_log); exit(0);}
	/**/
	/* Cells Num */
	cellnum=(xnumx-1)*(ynumy-1);
	if(cellnum>MAXCEL) {fprintf(fp_log,"Space out in pr[]"); fflush(fp_log); exit(0);}
	/**/
	/* Mark num */
	if(marknum>MAXMRK+1) {fprintf(fp_log,"Space out in markx[]"); fflush(fp_log); exit(0);}
	/**/
	/* Rock types Num */
	if (rocknum>MAXTMR){fprintf(fp_log,"Space out in marknu[]"); fflush(fp_log); exit(0);}
	/**/
	/* Bondary condit Equations Num */
	if (bondnum>MAXBON){fprintf(fp_log,"Space out in bondv[]"); fflush(fp_log); exit(0);}
	/**/
	/* Koef for processing */
	xstpx=xsize/(double)(xnumx-1);
	ystpy=ysize/(double)(ynumy-1);
	kfx=1.0/xstpx;
	kfy=1.0/ystpy;
	kfxx=kfx*kfx;
	kfyy=kfy*kfy;
	kfxy=kfx*kfy;
	/* Marker size */
	mardx=xstpx/(double)(mnumx);
	mardy=ystpy/(double)(mnumy);
	/* Step for numerical differentiation */
	numdx=5e-1*mardx;
	numdy=5e-1*mardy;
	/**/
	/* Spec counters */
	nodenum2=nodenum*2;
	nodenum3=nodenum*3;
	xnumx1=xnumx-1;
	ynumy1=ynumy-1;
	/**/
}


// Check if boundaries arrays are exceeded
void check_bound(int min,int max, int actual, const char error[])
{
	if (actual < min) { fprintf(fp_log,"ERROR: Array boundary check: %d < %d in %s \n",actual,min,error); fflush(fp_log); exit(0); }
	if (actual > max) { fprintf(fp_log,"ERROR: Array boundary check: %d > %d in %s \n",actual,max,error); fflush(fp_log); exit(0); }
}
