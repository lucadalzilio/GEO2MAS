/* Solve XY-Stokes+Continuity equations by vX,vY,P mixed arbitrary order Finite Difference method */
void viterateomp(int m0)
{
	/* Counters */
	long int m1,m2,m3,m1min,m1max,m2min,m2max,mcmax,mcmax0,mcmax1,dm1,dm2,ccc=0;
	double vxmin,vxmax,vymin,vymax,minvx,maxvx,minvy,maxvy,mindx,pmpa,minnu,maxnu,minpr,maxpr,mukoef,ridkoef,vel_cc,vel_cc0,vel_cc1,vel_cc2,timecc0,timecc1;
	double e,maxds,dss1,dsxx,dsxx1,dspp,dspp1,dsxy,dsxy1,dexx,dexx1,dexy,dexy1,celdx,celdy,swt1,dvx,dvy,dvx1,dvy1,mgg,mnu,mro,xelvis0=1.0,xelvis1=0;
	double dol0[MAXNOD],dolSxx[MAXNOD],dolSpp[MAXNOD],dolSxy[MAXNOD],dolExx[MAXNOD],dolExy[MAXNOD],dolVx[MAXNOD],dolVy[MAXNOD];
	int printyn=printmod,n1;
	/* Val buffer */
	double ival,res_low;
	int start_nx;
	/* Err koef */
	double bondsum,bondnum;
	double stoksum,stoknum;
	double contsum,contnum;
	double timestepe0;
	// Parallelization
	double start;
	int nt;
	long int m10,m20;
	double TK,EXY,EXYE,SXY,SXYE,EXX,SXX,PR,SXXE,SPPE,EXXE,VX,MVX,VY,MVY;

	start=omp_get_wtime();



	/* Change Grid, if e.g., need to re-arrange high resolution area */
	// Do only once; in mode.t3c set nr timestep/output to 1 to make prn, and than switch this off again
	int change_grid = 0;
#if setup > 9
	if(change_grid)
	{
		/* Save old erosion/sedimentation, hydration surfaces */
		for (m1=0;m1<xnumx;m1++)
		{
			ep0[m1]=ep[m1];
			ep0[xnumx+m1]=ep[xnumx+m1];
			ep0[xnumx*2+m1]=ep[xnumx*2+m1];
			ep0[xnumx*3+m1]=ep[xnumx*3+m1];
		}
	
		/* Recomputing horizontal grid */
		// shift_km in load input
		res_low = 2e3;
		start_nx = (int)(shift_km/res_low);
	
		ep0[xnumx*4] = 0;
		//  Shift all grid spacings to left, so HR area covers region of interest again
		//  C starts at 0, so gx[xnumx]=0 already
		for (m1=1;m1<xnumx-start_nx-1;m1++)
		{
			ep0[xnumx*4+m1]=ep0[xnumx*4+m1-1]+(gx[m1+start_nx+1]-gx[m1+start_nx]);
		}
		// Fill up gap at end with low resolution cells
		for (m1=xnumx-start_nx-1;m1<xnumx;m1++)
		{
			ep0[xnumx*4+m1]=ep0[xnumx*4+m1-1]+res_low;
		}

		//fprintf(fp_log,"%e %e %e %e %e %e %e %e", ep0[xnumx*4+0],ep0[xnumx*4+1],ep0[xnumx*4+251],ep0[xnumx*4+400],ep0[xnumx*4+401],ep0[xnumx*4+252],ep0[xnumx*4+xnumx-2],ep0[xnumx*4+xnumx-1]);getchar();
	
		/* Reinterpolate erosion/sedimentation, hydration surfaces */
		for (m3=0;m3<xnumx;m3++)
		{
			m1=m1serch(ep0[xnumx*4+m3]);
			/* Relative normalized coordinate calc */
			e=(ep0[xnumx*4+m3]-gx[m1])/(gx[m1+1]-gx[m1]);
		
			/* Surface level elevation for marker definition */
			ep[m3]=(e*ep0[m1+1]+(1.0-e)*ep0[m1]);
			ep[xnumx+m3]=(e*ep0[xnumx+m1+1]+(1.0-e)*ep0[xnumx+m1]);
			ep[xnumx*2+m3]=(e*ep0[xnumx*2+m1+1]+(1.0-e)*ep0[xnumx*2+m1]);
			ep[xnumx*3+m3]=(e*ep0[xnumx*3+m1+1]+(1.0-e)*ep0[xnumx*3+m1]);
		}
	
		/* Reload new gridline positions */
		for (m1=0;m1<xnumx;m1++)
		{
			gx[m1]=ep0[xnumx*4+m1];
		}
	
		/* Reload High Resolution area limits */
		// m10_hr - m11_hr
		// m20_hr
		//   |
		// m21_hr		
		res_high = 100.0e3;
		for (m1=0;m1<xnumx;m1++)
		{
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
	
	}
#endif
	/* End Change Grid */


	/* -- Step 2: Calculation of scalar and tensor properties: preparation for solving mechanical eq. -- */

	/* Recomputing marker properties with new timestep */
	if (stoksmod)
	{
		if (printmod) fprintf(fp_log,"\n RO, NU, CP etc  RECALC..."); fflush(fp_log);
		
		// ! Only effective call to ronurecalc! Here will go to viscalc !
		ronurecalcomp();
		
		if (printmod) fprintf(fp_log," RO, NU, CP etc OK!\n"); fflush(fp_log);
		if (printmod) fprintf(fp_log,"Elastic timestep %e \n",timestepe/3.15576e+7); fflush(fp_log);
	}

	/* -- Step 3: Solving momentum and continuity equations -- */

	/* Step 3.1. Check timesteps and print values   */
	if(stoksmod)
	{
		/* Movement Timestep check   */
		if (timestep<=0)
		{
			fprintf(fp_log,"EXIT PROGRAM:  Movement timestep<=0 <%e>",timestep);fflush(fp_log);
			exit(0);
		}
		/* Elastic Timestep check   */
		if (timestepe<=0)
		{
			fprintf(fp_log,"EXIT PROGRAM:  Elastic timestep<=0 <%e>",timestepe);fflush(fp_log);
			exit(0);
		}

		/* Check  Xelvis by printing the min and max of this Z: find min and max */	
		/* Temporarily change stoksmod for Xelvis calc mode, so can print to log file */	
		stoksmod=-stoksmod;
		for (m1=1;m1<xnumx;m1++)
		{
			for (m2=1;m2<ynumy;m2++)
			{
				/* Pos in vx[], vy[], pr[], etc. */
				mcmax1=m1*ynumy+m2;

				/* Sxx */	
				ival=sxxcalc(m1,m2,0);
				if(ival<xelvis0) xelvis0=ival; if(ival>xelvis1) xelvis1=ival;

				/* Sxy,Exy */	
				if(m1<xnumx-1 && m2<ynumy-1)
				{
					ival=sxycalc(m1,m2,0);
					if(ival<xelvis0) xelvis0=ival; if(ival>xelvis1) xelvis1=ival;
				}
			}
		}
		
		/* Restore stoksmod */	
		stoksmod=-stoksmod;
		if(printmod) fprintf(fp_log,"\n !!!   VISCOELASTIC TIME STEP FOR CYCLE %e YEAR !!!\n",timestepe/3.15576e+7);
		if(printmod) fprintf(fp_log,"\n !!!   XELVIS:    min=%e    max=%e !!!\n",xelvis0,xelvis1);fflush(fp_log);
	}
	/* End Defining effective numerical timestep for Stockes Equation */

	/* Step 3.2. Clear New Solution */
	/* P, Vx,Vy */
	pos0cur=0;
	/* Err koef */
	bondsum=0;bondnum=0;
	stoksum=0;stoknum=0;
	contsum=0;contnum=0;
	timestepe0=timestepe;

	/* Step 3.3. Time step for markers definition for drunken sailor instability correction */
	if(markmod)
	{
		// Vx,Vy max-min definition 
		minvx=1e+30;maxvx=-1e+30;minvy=1e+30;maxvy=-1e+30;
		for (m1=0;m1<xnumx;m1++)
			for (m2=0;m2<ynumy;m2++)
		{
			// Pos of Vx in sol0[] 
			m3=m1*ynumy+m2;
                
			// Min,Max Vx definition 
			if(m2<ynumy-1)
			{
				minvx=MINV(minvx,vx[m3]);
				maxvx=MAXV(maxvx,vx[m3]);
			}
			// Min,Max Vy definition 
			if(m1<xnumx-1)
			{
				minvy=MINV(minvy,vy[m3]);
				maxvy=MAXV(maxvy,vy[m3]);
			}
		}
		// Max Vx,Vy Diff in Grid Calc 
		vxmin=minvx;
		vymin=minvy;
		vxmax=maxvx;
		vymax=maxvy;
		maxvx=MAXV(ABSV(vxmin),ABSV(vxmax));
		maxvy=MAXV(ABSV(vymin),ABSV(vymax));
		if (maxvx)
		{
			maxvx=(maxxystep*res_high)/maxvx;
			if(printmod) fprintf(fp_log,"\n !!! MAX VALID TIME STEP FOR Vx-MARKER %e YEAR !!!\n",maxvx/3.15576e+7); fflush(fp_log);
			if(timestep>maxvx) timestep=maxvx;
		}
		if (maxvy)
		{
			maxvy=(maxxystep*res_high)/maxvy;
			if(printmod) fprintf(fp_log,"\n !!! MAX VALID TIME STEP FOR Vy-MARKER %e YEAR !!!\n",maxvy/3.15576e+7); fflush(fp_log);
			if(timestep>maxvy) timestep=maxvy;
		}
		if(printmod) fprintf(fp_log,"\n !!!       CURRENT INITIAL DISPLACEMENT TIME STEP FOR CYCLE %e YEAR !!!\n",timestep/3.15576e+7); fflush(fp_log);
	}

	/* Step 3.4. Add Matrix by vX-vY-Stokes, Continuity, Boundary, EPS, SIG, Equations */
	mukoef=nubeg*strmin;
	/* Node  Cycle */
	for (m1=0;m1<xnumx;m1++)
		for (m2=0;m2<ynumy;m2++)
	{
		/* Pos P,Vx,Vy in sol0[] */
		mcmax=(m1*ynumy+m2)*3;

		/* Add Continuity equation for Cells ------------------------------------------------ */
		if(m1 && m2)
		{
			if(!bondm[mcmax+0]) 
			{
				/* Continuity eq. */
				conterr(m1,m2,0);
			}
			else
			{
				/* Add P-Boundary */
				pbonderr(mcmax+0,0);
			}
			/* Rescale coefficients */
			for (n1=0;n1<=wn[0];n1++) wi[n1]*=mukoef;
			gausmat4(1,mcmax+0,0);
		}
		else
		{
			wn[0]=1;
			wi[0]=0;
			wn[1]=mcmax+0;
			wi[1]=mukoef;
			gausmat4(1,mcmax+0,0);
		}
		/* Add vX-Equations --------------------------------------------------- */
		if(m2<ynumy-1)
		{
			if(!bondm[mcmax+1] || (timesum>timebond && m1>2 && m2>2 && m1<xnumx-4 && m2<ynumy-3)) 
			{
				/* Add vX-Stokes */
				// One could hardcode a runtime change of push velocity here 
				// Remember: dx=2km, dy=0.5km
			
				xstokserr(m1,m2,0);
				/**/
				/* Add matrix: vX */
				gausmat4(1,mcmax+1,0);
			}
			else
			{
				/* Continuity Equation Vx boundary condition */
				if(bondv[bondm[mcmax+1]][1]>1e+30)
				{
					/* m1 m2 increment definition */
					m3=(bondn[bondm[mcmax+1]][0]-1-(mcmax+1))/3;
					dm1=dm2=0;
					if(m3>=ynumy) dm1=1;
					if(m3==1 || m3>ynumy) dm2=1;
					if(m3<0 || m3>ynumy+1 || !(m1+dm1) || !(m2+dm2) || !bondm[bondn[bondm[mcmax+1]][0]-2])
					{
						fprintf(fp_log,"EXIT PROGRAM Inconsistent Vx(%ld) P(%ld) boundary conditions for X=%ld Y=%ld Node",mcmax+1,bondn[bondm[mcmax+1]][0]-2,m1,m2);
						fflush(fp_log);
						exit(0);
					}
					conterr(m1+dm1,m2+dm2,0);
				}
				else
				{
					/* Add vX Simple Boundary */
					xbonderr(mcmax+1,0);
				}
				
				/* Vx boundary condition Add */
				/* Rescale coefficients */
				for (n1=0;n1<=wn[0];n1++) wi[n1]*=mukoef;
				gausmat4(1,mcmax+1,0);
			}
		}
		/* Add Ghost parameters */
		else
		{
			wn[0]=1;
			wi[0]=0;
			wn[1]=mcmax+1;
			wi[1]=mukoef;
			gausmat4(1,mcmax+1,0);
		}

		/* Add vY-Equations --------------------------------------------------- */
		if(m1<xnumx-1)
		{
			if(!bondm[mcmax+2]) 
			{
				/* Add vX-Stokes */
				ystokserr(m1,m2,0);

				/* Add matrix: vY */
				gausmat4(1,mcmax+2,0);
			}
			else
			{
				/* Continuity Equation Vy boundary condition */
				if(bondv[bondm[mcmax+2]][1]>1e+30)
				{
					/* m1 m2 increment definition */
					m3=(bondn[bondm[mcmax+2]][0]-1-(mcmax+2))/3;
					dm1=dm2=0;
					if(m3>=ynumy) dm1=1;
					if(m3==1 || m3>ynumy) dm2=1;
					if(m3<0 || m3>ynumy+1 || !(m1+dm1) || !(m2+dm2) || !bondm[bondn[bondm[mcmax+2]][0]-3])
					{
						fprintf(fp_log,"EXIT PROGRAM Inconsistent Vy(%ld) P(%ld) boundary conditions for X=%ld Y=%ld Node",mcmax+2,bondn[bondm[mcmax+2]][0]-3,m1,m2);
						fflush(fp_log);
						exit(0);
					}
					conterr(m1+dm1,m2+dm2,0);
				}
				/* Simple Vy boundary condition */
				else
				{
					/* Add vY Simple Boundary */
					ybonderr(mcmax+2,0);
				}

				/* Vy boundary condition Add */
				/* Rescale coefficients */
				for (n1=0;n1<=wn[0];n1++) wi[n1]*=mukoef;
				gausmat4(1,mcmax+2,0);
			}
		}
		/* Add Ghost parameters */
		else
		{
			wn[0]=1;
			wi[0]=0;
			wn[1]=mcmax+2;
			wi[1]=mukoef;
			gausmat4(1,mcmax+2,0);
		}
	}
	/* End  Add Matrix By vX-vY-Stokes, Continuity Equations */

	if (printmod==10000) fprintf(fp_log," Time taken for adding coefficients of Stokes+Continuity equations = %e s \n",omp_get_wtime()-start);
	
	/* Solve Matrix ================================================ */
	start=omp_get_wtime();
	if (printmod) fprintf(fp_log,"Number of positions in global matrix = %ld  Number of unknown = %ld \n",pos0cur,nodenum3); fflush(fp_log);
	gausmat4(0,nodenum3,0);
	if (printmod==10000) fprintf(fp_log," Time taken for solving Stokes+Continuity equations = %e s ",omp_get_wtime()-start);
	/* Solve Matrix ================================================ */

	start=omp_get_wtime();
	/* Reload P, Vx, Vy Results */
	/* Node  Cycle */
	for (m1=0;m1<xnumx;m1++)
	{
		/* Set Initial p value at upper boundary */
		pmpa=pinit;

		for (m2=0;m2<ynumy;m2++)
		{
			/* Pos P,Vx,Vy in sol0[] */
			mcmax0=(m1*ynumy+m2)*3;
			/* Pos in vx[], vy[], pr[], etc. */
			mcmax1=m1*ynumy+m2;

			/* Reload/Recalc P */	
			if(m1 && m2) 
			{
				/* Reload P */	
				pr[mcmax1]=x[mcmax0+0];
				/*	
				fprintf(fp_log,"\n %ld %e",mcmax1,pr[mcmax1]); getchar();
				*/	
			}

			/* Reload Vx */	
			if(m2<ynumy-1)
			{ 
				vx[mcmax1]=x[mcmax0+1];
			}
	
			/* Reload Vy */	
			if(m1<xnumx-1)
			{
				vy[mcmax1]=x[mcmax0+2];
			}
		}
	}
	/* End Reload P, Vx, Vy Results */
	
	/* Recalc EPS, SIG Results */
	/* Node  Cycle */
	maxds=0;
	minnu=1e+50;maxnu=-1e+50;
	for (m1=1;m1<xnumx;m1++)
		for (m2=1;m2<ynumy;m2++)
	{
		/* Pos in vx[], vy[], pr[], etc. */
		mcmax1=m1*ynumy+m2;
	
		/* Sxx,Exx */	
		sxx[mcmax1]=sxxcalc(m1,m2,0); exx[mcmax1]=eps[0];
		maxds=MAXV(maxds,ABSV(sxx[mcmax1]-sxx0[mcmax1]));
		/* Min,Max Nu value Calc */
		minnu=MINV(minnu,eps[2]);
		maxnu=MAXV(maxnu,eps[2]);

		/* Sxy,Exy */	
		if(m1<xnumx-1 && m2<ynumy-1)
		{
			sxy[mcmax1]=sxycalc(m1,m2,0); exy[mcmax1]=eps[0]; esp[mcmax1]=eps[1];
			maxds=MAXV(maxds,ABSV(sxy[mcmax1]-sxy0[mcmax1]));
		}
	}
	/* End Recalc EPS, SIG Results */

	/* Vx,Vy max-min definition */
	minvx=1e+30;maxvx=-1e+30;minvy=1e+30;maxvy=-1e+30;
	for (m1=0;m1<xnumx;m1++)
		for (m2=0;m2<ynumy;m2++)
	{
		/* Pos of Vx in sol0[] */
		m3=m1*ynumy+m2;
		/**/
		/* Min,Max Vx definition */
		if(m2<ynumy-1)
		{
			minvx=MINV(minvx,vx[m3]);
			maxvx=MAXV(maxvx,vx[m3]);
		}
		/* Min,Max Vy definition */
		if(m1<xnumx-1)
		{
			minvy=MINV(minvy,vy[m3]);
			maxvy=MAXV(maxvy,vy[m3]);
		}
	}
	/* Max Vx,Vy Diff in Grid Calc */
	vxmin=minvx;
	vymin=minvy;
	vxmax=maxvx;
	vymax=maxvy;
	maxvx-=minvx;
	maxvy-=minvy;
	maxvx=MAXV(maxvx,maxvy);
	mindx=res_high;

	/* Check Error */
	/* Node  Cycle */
	minpr=1e+50;maxpr=-1e+50;
	for (m1=0;m1<xnumx;m1++)
		for (m2=0;m2<ynumy;m2++)
	{
		/* Pos of P,Vx,Vy in sol0[] */
		mcmax0=(m1*ynumy+m2)*3;

		/* Check Continuity equation for Cells =========================== */
		if(m1 && m2) 
		{
			ival=conterr(m1,m2,1);
			/*	
			fprintf(fp_log,"\n %ld %ld   %e %e %e     %e ",m1,m2,ival,maxvx,mindx,ival/(maxvx/mindx)); getchar();
			*/	
			contsum+=ival*ival;
			contnum+=1.0;
			/* Print Results */
			ival/=maxvx/mindx;
			if (printmod && ABSV(ival)>DIVVMIN)
			{
				fprintf(fp_log,"\n Large Continuity err at X=%ld Y=%ld:   Err=%e",m1,m2,ival); fflush(fp_log);
			}
		}

		/* Check vX-Equations for nodes =========================== */
		if(m2<ynumy-1)
		{
			if(!bondm[mcmax0+1]) 
			{
				/* Add vX-Stokes */
				ival=xstokserr(m1,m2,1);
				stoksum+=ival*ival;
				stoknum+=1.0;
				/* Min,Max Pr value Calc */
				maxpr=MAXV(maxpr,errbuf[10]);
				minpr=MINV(minpr,errbuf[11]);
				/* Print Results */
				ival/=maxvx/mindx/mindx;
				if (printmod && ABSV(ival)>STOKSMIN)
				{
					fprintf(fp_log,"\n Large X stokes err at X=%ld Y=%ld:   Err=%e",m1,m2,ival); fflush(fp_log);
				}
			}
			else
			{
				if(bondv[bondm[mcmax0+1]][1]<1e+30)
				{
					/* Add vX-Boundary */
					ival=xbonderr(mcmax0+1,1);
					bondsum+=ival*ival;
					bondnum+=1.0;
				}
			}
		}

		/* Check vY-Equations for nodes =========================== */
		if(m1<xnumx-1)
		{
			if(!bondm[mcmax0+2]) 
			{
				/* Add vX-vY-Stokes */
				ival=ystokserr(m1,m2,1);
				stoksum+=ival*ival;
				stoknum+=1.0;
				/* Min,Max Pr value Calc */
				maxpr=MAXV(maxpr,errbuf[10]);
				minpr=MINV(minpr,errbuf[11]);
				/* Print Results */
				ival/=maxvx/mindx/mindx;
				if (printmod && ABSV(ival)>STOKSMIN)
				{
					fprintf(fp_log,"\n Large Y stokes err at X=%ld Y=%ld:   Err=%e",m1,m2,ival); fflush(fp_log);
				}
			}
			else
			{
				if(bondv[bondm[mcmax0+2]][1]<1e+30)
				{
					/* Add vX-vY-Boundary */
					ival=ybonderr(mcmax0+2,1);
					bondsum+=ival*ival;
					bondnum+=1.0;
				}
			}
		}
	}
	stoksum=pow(stoksum/stoknum,0.5)/(maxvx/mindx/mindx);
	contsum=pow(contsum/contnum,0.5)/(maxvx/mindx);
	bondsum=pow(bondsum/bondnum,0.5)/maxvx;
	/* End Check Error */

	/* Print Results */
	if (printmod)
	{
		fprintf(fp_log,"\n KRUG %2d \n MIN/MAX VELOCITY Vx = % e / % e     Vy = %e / % e\n",m0+1,vxmin,vxmax,vymin,vymax);
		fprintf(fp_log,"VISKOS: min = %e max = %e \n",minnu,maxnu);
		fprintf(fp_log,"PRESSURE: min = %e max = %e \n",minpr,maxpr);
		fprintf(fp_log,"STRESS CHANGE:    max = %e \n",maxds);
		fprintf(fp_log,"STOKS: num = %e err = %e \n",stoknum,stoksum);
		fprintf(fp_log,"CONT : num = %e err = %e \n",contnum,contsum);
		fprintf(fp_log,"BOUND V: num = %e err = %e \n",bondnum,bondsum);
		fflush(fp_log);
	}



	/* --> Step 4: Set time step for marker displacement <-- */
	
	timestep=maxtmstep;
	
	// Check if markers would not be displaced beyond set limit of maxxystep set in mode.t3c
	if(markmod)
	{
		maxvx=MAXV(ABSV(vxmin),ABSV(vxmax));
		maxvy=MAXV(ABSV(vymin),ABSV(vymax));
		if (maxvx)
		{
			maxvx=(maxxystep*res_high)/maxvx;
			if(printmod) fprintf(fp_log,"\n !!! MAX VALID TIME STEP FOR Vx-MARKER %e YEAR !!!\n",maxvx/3.15576e+7); fflush(fp_log);
			if(timestep>maxvx) timestep=maxvx;
		}
		if (maxvy)
		{
			maxvy=(maxxystep*res_high)/maxvy;
			if(printmod) fprintf(fp_log,"\n !!! MAX VALID TIME STEP FOR Vy-MARKER %e YEAR !!!\n",maxvy/3.15576e+7); fflush(fp_log);
			if(timestep>maxvy) timestep=maxvy;
		}
		if(printmod) fprintf(fp_log,"\n !!!       CURRENT TIME STEP FOR CYCLE %e YEAR !!!\n",timestep/3.15576e+7); fflush(fp_log);
	}

	/* Recalc EPS, SIG Results for new smaller displacement timestep  */
	/* Node  Cycle */
	if(stoksmod && timestep<timestepe)
	{
		// Set timestepe to displacement time step instead, so when into sxxcalc will get stresses for new displacement time step following same non-linear equation. Test wether simple linear interpolation with new stress with computational time step is better.
		timestepe=timestep;
		
		maxds=0;
		if(printmod) fprintf(fp_log,"\n !!! RECALC SIGMA FOR NEW TIMESTEP=%e YR !!!\n",timestep/3.15576e+7); fflush(fp_log);
		for (m1=1;m1<xnumx;m1++)
			for (m2=1;m2<ynumy;m2++)
		{
			/* Pos in vx[], vy[], pr[], etc. */
			mcmax1=m1*ynumy+m2;
			/**/	
			/* Sxx,Exx */	
			sxx[mcmax1]=sxxcalc(m1,m2,0); exx[mcmax1]=eps[0];
			maxds=MAXV(maxds,ABSV(sxx[mcmax1]-sxx0[mcmax1]));
			/**/	
			/* Sxy,Exy */	
			if(m1<xnumx-1 && m2<ynumy-1)
			{
				sxy[mcmax1]=sxycalc(m1,m2,0); exy[mcmax1]=eps[0]; esp[mcmax1]=eps[1];
				maxds=MAXV(maxds,ABSV(sxy[mcmax1]-sxy0[mcmax1]));
			}
		}
		if(printmod) fprintf(fp_log,"NEW STRESSS CHANGE:    max = %e \n",maxds); fflush(fp_log);
	}
	/* End Recalc EPS, SIG Results for new timestep  */



	/* --> Step 5: Interpolate stress changes etc from nodes to markers <-- */ 
	
	/* Set values for NEWLY coming marker, i.e. T=0 or near BC where need to relax perturbations */
	for (int m7=0;m7<=marknum;m7++) 
	{
		/* If now come inside, NEW = had no temperature */
		if (markx[m7]>0 && marky[m7]>0 && markx[m7]<xsize && marky[m7]<ysize && markt[m7]<50 && (markk[m7]<=0 || (markx[m7]>200e3 && markx[m7]<450e3 && marky[m7]<85e3) ))
			//if (markx[m7]>0 && marky[m7]>0 && markx[m7]<xsize && marky[m7]<ysize && markt[m7]<50 && (markk[m7]<=0 || (markx[m7]>200e3 && markx[m7]<450e3 && marky[m7]<85e3) || (markx[m7]>1250e3 && markx[m7]<1440e3 && marky[m7]<85e3)))
		{
			/* Interpolate Stresses */
			m10=m1serch(markx[m7]);
			m20=m2serch(marky[m7]);
			allinterdomp(markx[m7],marky[m7],m10,m20,&TK,&EXY,&EXYE,&SXY,&SXYE,&EXX,&SXX,&PR,&SXXE,&SPPE,&EXXE,&VX,&MVX,&VY,&MVY);
			/**/
			/* For newly coming markers set old values of marker variables */
			markxx[m7]=SXXE;	// E = Old stresses, strain rates, and pressure
			markxy[m7]=SXYE;
			markp[m7]=SPPE;
			markexx[m7]=EXXE;
			markexy[m7]=EXYE;
			markvx[m7]=MVX;
			markvy[m7]=MVY;
		
			/* Set newly coming marker if no temperature calculation */
			if(!tempmod) markk[m7]=TK;
			/*
			markexx[m7]=EXX;
			markexy[m7]=EXY;
			fprintf(fp_log,"DIFFUSION1 %ld %e",m7,markk[m7]); getchar();
			*/
		}
	}

	/* Recalc marker SIGij + Diffusion, Add  nodes wt */
	
	// Prepare arrays for shared memory storage
	/* Clear nodes wt, save increment */
	for (m1=0;m1<nodenum;m1++)
	{
		dol0[m1]=0;
		dolSxx[m1]=0;
		dolSpp[m1]=0;
		dolSxy[m1]=0;
		dolExx[m1]=0;
		dolExy[m1]=0;
		dolVx[m1]=0;
		dolVy[m1]=0;
	}
	
#pragma omp parallel shared(markx,marky,markt,markk,markv,markxx,markp,markxy,markexx,markexy,markvx,markvy,markd,timestep,timestepe0,stredif,marknum,xsize,ysize,gx,gy,xnumx,ynumy,markgg,dol0,dolSxx,dolSpp,dolSxy,dolExx,dolExy,dolVx,dolVy) \
	private(m10,m20,ival,TK,SXX,PR,SXY,EXX,EXY,SXXE,SPPE,SXYE,EXXE,EXYE,VX,MVX,VY,MVY,dsxx,dspp,dsxy,dexx,dexy,dvx,dvy,dsxx1,dspp1,dsxy1,dexx1,dexy1,dvx1,dvy1,celdx,celdy,swt1,mgg,mnu,mro,dss1) 
	{
		
		// For each core keep a private array to be added latter under critical, since reduction clause does not work on arrays
		double* private_dol0  ;
		double* private_dolSxx;
		double* private_dolSpp;
		double* private_dolSxy;
		double* private_dolExx;
		double* private_dolExy;
		double* private_dolVx ;
		double* private_dolVy ;
		if (stredif) 
		{
			private_dol0   = (double*) calloc(nodenum,sizeof(double));	
			private_dolSxx = (double*) calloc(nodenum,sizeof(double));	
			private_dolSpp = (double*) calloc(nodenum,sizeof(double));	
			private_dolSxy = (double*) calloc(nodenum,sizeof(double));	
			private_dolExx = (double*) calloc(nodenum,sizeof(double));	
			private_dolExy = (double*) calloc(nodenum,sizeof(double));	
			private_dolVx  = (double*) calloc(nodenum,sizeof(double));	
			private_dolVy  = (double*) calloc(nodenum,sizeof(double));	
		}

#pragma omp for \
		schedule(runtime)
		// For all markers that were already inside
		for (long int m6=0;m6<=marknum;m6++) 
		{
			/* Check that exclude non-rock markers markers and those out of grid */
			if (markx[m6]>0 && marky[m6]>0 && markx[m6]<xsize && (marky[m6])<ysize && markt[m6]<50 && markk[m6]>0)
			{
				/* Interpolate Stresses */
				m10=m1serch(markx[m6]);
				m20=m2serch(marky[m6]);
				allinterdomp(markx[m6],marky[m6],m10,m20,&TK,&EXY,&EXYE,&SXY,&SXYE,&EXX,&SXX,&PR,&SXXE,&SPPE,&EXXE,&VX,&MVX,&VY,&MVY);

				/* Numerical diffusion add to markers -------------*/
				if(stredif)
				{
					/* Calc difference in old marker stresses, wether on nodes or markers  */
					dsxx=SXXE-markxx[m6];
					dspp=SPPE-markp[m6];
					dsxy=SXYE-markxy[m6];
					dexx=EXXE-markexx[m6];
					dexy=EXYE-markexy[m6];
					dvx=MVX-markvx[m6];
					dvy=MVY-markvy[m6];

					/* Interpolation of G and Nu from nodes to marker */
					/* Marker weight calculation using dimension of current Cell */
					celdx=gx[m10+1]-gx[m10];
					celdy=gy[m20+1]-gy[m20];
					swt1=1.0/celdx/celdy;

					/* Calc, check Viscous relaxation to marker stresses */
					mgg=markgg[markt[m6]];
					mnu=markv[m6];
					mro=markd[m6];
			
					// Calculate diffusion term	
					/* dS=dS0*exp(-G/Nu*dt) */
					dss1=-stredif*mgg/mnu*timestep;
					if(dss1<-150.0) dss1=-150.0;
					dss1=1.0-exp(dss1);
					dsxx1=dsxx*dss1;
					dspp1=dspp*dss1;
					dsxy1=dsxy*dss1;
					dexx1=dexx*dss1;
					dexy1=dexy*dss1;
		
					/* dV=dV0*exp(-Nu/Ro*dt*(2/dx^2+2/dy^2)) */
					dss1=-stredif*mnu/mro*(2.0/celdx/celdx+2.0/celdy/celdy)*timestep;
					if(dss1<-150.0) dss1=-150.0;
					dss1=1.0-exp(dss1);
					dvx1=dvx*dss1;
					dvy1=dvy*dss1;

					/* Diffuse Marker properties */
					markxx[m6]+=(dsxx1);
					markp[m6]+=(dspp1);
					markxy[m6]+=(dsxy1);
					markexx[m6]+=(dexx1);
					markexy[m6]+=(dexy1);
					markvx[m6]+=(dvx1);
					markvy[m6]+=(dvy1);

					/* Wt for nodes calc, add */
					celdx=(markx[m6]-gx[m10])/(gx[m10+1]-gx[m10]);
					celdy=(marky[m6]-gy[m20])/(gy[m20+1]-gy[m20]);
					if (celdx<0 || celdx>1 || celdy <0 || celdy>1 || swt1<0 ) {fprintf(fp_log,"ERROR @ weigths in diffusion: %ld %e %e %e \n",m6,celdx,celdy,swt1); fflush(fp_log); getchar();}
			
					for (int m1=m10;m1<=m10+1;m1++)
					{	
						for (int m2=m20;m2<=m20+1;m2++)
						{
							/* Cur Node Num, wt */
							long int m4=m1*ynumy+m2;
							if (m1==m10 && m2==m20)
							{
								ival=(1.0-celdx)*(1.0-celdy)*swt1;
							}
							else if (m1==m10+1 && m2==m20)
							{
								ival=(celdx)*(1.0-celdy)*swt1;
							}
							else if (m1==m10 && m2==m20+1)
							{
								ival=(1.0-celdx)*(celdy)*swt1;
							}
							else if (m1==m10+1 && m2==m20+1)
							{
								ival=(celdx)*(celdy)*swt1;
							}	

							/* Add Node wt, T */
							private_dol0[m4]+=ival;
							private_dolSxx[m4]+=dsxx1*ival;
							private_dolSpp[m4]+=dspp1*ival;
							private_dolSxy[m4]+=dsxy1*ival;
							private_dolExx[m4]+=dexx1*ival;
							private_dolExy[m4]+=dexy1*ival;
							private_dolVx[m4]+=dvx1*ival;
							private_dolVy[m4]+=dvy1*ival;
						}
					}
				}
				/* End - Numerical diffusion added to markers -------------*/

				/* Change marker stresses after solution */
				markxx[m6]+=(SXX-SXXE);
				markp[m6]+=(PR-SPPE);
				markxy[m6]+=(SXY-SXYE);
				markexx[m6]+=(EXX-EXXE);
				markexy[m6]+=(EXY-EXYE);
				markvx[m6]+=(VX-MVX)*timestep/timestepe0;
				markvy[m6]+=(VY-MVY)*timestep/timestepe0;
			}
		}

		// Summation of arrays from different processors
		if(stredif)
		{
#pragma omp critical (sumdolarrays)
			{
				for (long int m8=0;m8<=nodenum;m8++)
				{
					dol0[m8]+=private_dol0[m8];
					dolSxx[m8]+=private_dolSxx[m8];
					dolSpp[m8]+=private_dolSpp[m8];
					dolSxy[m8]+=private_dolSxy[m8];
					dolExx[m8]+=private_dolExx[m8];
					dolExy[m8]+=private_dolExy[m8];
					dolVx[m8]+=private_dolVx[m8];
					dolVy[m8]+=private_dolVy[m8];
				}
			}
	
			// Free allocated dynamic dol arrays
			free(private_dol0);
			free(private_dolSxx);
			free(private_dolSpp);
			free(private_dolSxy);
			free(private_dolExx);
			free(private_dolExy);
			free(private_dolVx);
			free(private_dolVy);
		}
	}
	// End OMP section
	
	/* Numerical antidiffusion add to markers -------------*/
	if(stredif)
	{
		/* Recalc changes in nodes T */
		for (m1=0;m1<nodenum;m1++)
		{
			if(dol0[m1]) 
			{
				/* Averaged Changes in Stresses due to smoothing */
				dolSxx[m1]/=dol0[m1];
				dolSpp[m1]/=dol0[m1];
				dolSxy[m1]/=dol0[m1];
				dolExx[m1]/=dol0[m1];
				dolExy[m1]/=dol0[m1];
				dolVx[m1]/=dol0[m1];
				dolVy[m1]/=dol0[m1];
			}
		}

		/* Recalc marker SIGij + Antidiffusion) */
		// shared by default now; so if need have; copyin(tk,exy,exye,exx,exxe,sxy,sxye,sxx,sxxe,pr,sppe,vx,mvx,vy,mvy) 
#pragma omp parallel for shared(markx,marky,markt,markk,markxx,markp,markxy,markexx,markexy,markvx,markvy,dolSxx,dolSpp,dolSxy,dolExx,dolExy,dolVx,dolVy,marknum,xsize,ysize,gx,gy,xnumx,ynumy) \
		private(m10,m20,ival,celdx,celdy) \
		schedule(runtime)
		for (long int m6=0;m6<=marknum;m6++) 
		{
			/* Check that exclude non-rock markers markers and those out of grid */
			if (markx[m6]>0 && marky[m6]>0 && (markx[m6])<xsize && (marky[m6])<ysize && markt[m6]<50 && markk[m6]>0)
			{
				
				m10=m1serch(markx[m6]);
				m20=m2serch(marky[m6]);
				
				/* Wt for nodes calc, add */
				celdx=(markx[m6]-gx[m10])/(gx[m10+1]-gx[m10]);
				celdy=(marky[m6]-gy[m20])/(gy[m20+1]-gy[m20]);
				if (celdx<0 || celdx>1 || celdy <0 || celdy>1) {fprintf(fp_log,"ERROR @ weights antidiffusion: %ld %e %e %e \n",m6,celdx,celdy,swt1); fflush(fp_log); getchar();}
		
				for (int m1=m10;m1<=m10+1;m1++)
				{
					for (int m2=m20;m2<=m20+1;m2++)
					{
						/* Cur Node Num, wt */
						long int m4=m1*ynumy+m2;
			
						// calculate weight from distances
						if (m1==m10 && m2==m20)
							{ival=(1.0-celdx)*(1.0-celdy);}
						else if (m1==m10+1 && m2==m20)
							{ival=(celdx)*(1.0-celdy);}
						else if (m1==m10 && m2==m20+1)
							{ival=(1.0-celdx)*(celdy);}
						else if (m1==m10+1 && m2==m20+1)
							{ival=(celdx)*(celdy);}
       		 	
						/* Antidiffuse Marker Stresses */
						if(dolSxx[m4]) markxx[m6]-=dolSxx[m4]*ival;
						if(dolSpp[m4]) markp[m6]-=dolSpp[m4]*ival;
						if(dolSxy[m4]) markxy[m6]-=dolSxy[m4]*ival;
						if(dolExx[m4]) markexx[m6]-=dolExx[m4]*ival;
						if(dolExy[m4]) markexy[m6]-=dolExy[m4]*ival;
						if(dolVx[m4]) markvx[m6]-=dolVx[m4]*ival;
						if(dolVy[m4]) markvy[m6]-=dolVy[m4]*ival;
					}
				}
			}
		}
	}
	/* End Numerical antidiffusion add to markers -------------*/
	
	if (printmod==10000) fprintf(fp_log," Time taken final part viterate (viscoelastic stress calc, (anti)diffusion ) = %e s \n",omp_get_wtime()-start); fflush(fp_log);
}
/* End Solve XY-Stokes+Continuity equations by vX,vY,P mixed arbitrary order Finite Difference method */



/* LEFT+Right Side or Err for X-Stokes Equation */
/* Stoks equation initial form */
/* dSIGxx/dX + dSIGxy/dY - dP/dX = RO*DVx/Dt - RO*Gx */
double xstokserr(long int m1, long int m2, int ynerr)
	/* m1,m2 - node X,Y number */
	/* ynerr - Err Calc Y(1)/N(0) */
{
	/* Counters */
	int n1,n2;
	long int v[4];
	/* Err Buf */
	double leftx,rightx,nueff;
	/* Distances */
	double xkf=(gx[m1+1]-gx[m1-1])/2.0,ykf=gy[m2+1]-gy[m2];

	/* Staggered Nodes num */
	/*                      [0]                [2]   */
	/*                  Ro0,Sxy0,Nu0                 */
	/*                                               */
	/*         Pr1,Sxx1    <Vx0>  Pr3,Sxx3           */
	/*                                               */
	/*                      [1]                [3]   */
	/*                  RO1,Sxy1,Nu1                 */
	/*                                               */
	v[0]=m1*ynumy+m2;v[1]=v[0]+1;
	v[2]=v[0]+ynumy;v[3]=v[2]+1;

	/* RIGHT parts of X-Stokes */
	/* dSIGxx/dX + dSIGxy/dY - dP/dX = RO*DVx/Dt - RO*Gx */
	rightx  = (-GXKOEF-inertyn*mvx[v[0]]/timestepe)*mrx[v[0]];

	/* Return val for LSQ err ----------------------------*/
	if(ynerr==1)
	{
		/* LEFT part of X-Stokes */
		/* dSIGxx/dX + dSIGxy/dY - dP/dX = - RO*Gx */

		/* dSIGxx/dX */
		leftx =(sxx[v[3]]-sxx[v[1]])/xkf;
		/* dSIGxy/dY */
		leftx+=(sxy[v[1]]-sxy[v[0]])/ykf;
		/* -dP/dX */
		leftx-=(pr[v[3]]-pr[v[1]])/xkf;
		/* DVx/Dt*RO */
		leftx-=inertyn*vx[v[0]]*mrx[v[0]]/timestepe;

		/* Next 2 lines added for drunken sailor instability */
		/* -gx/dt/2*(vx*dRHO/dx+vy*dRHO/dy) */
		leftx-=timestep/2.0*GXKOEF*(vx[v[0]]*(ro[v[2]]+ro[v[3]]-ro[v[0]-ynumy]-ro[v[1]-ynumy])/(gx[m1+1]-gx[m1-1])/2.0+((vy[v[0]]+vy[v[1]])*(gx[m1]-gx[m1-1])+(vy[v[0]-ynumy]+vy[v[1]-ynumy])*(gx[m1+1]-gx[m1]))/(gx[m1+1]-gx[m1-1])/2.0*(ro[v[1]]-ro[v[0]])/(gy[m2+1]-gy[m2]));

		/* Effective NU calc */
		nueff=MAXV(ABSV(nu[v[0]]),ABSV(nu[v[1]]));

		/* X-STOKS Error */
		leftx=(leftx-rightx)/nueff;

		/* Min,Max Value of P Save */
		errbuf[10]=MAXV(pr[v[1]],pr[v[3]]);
		errbuf[11]=MINV(pr[v[1]],pr[v[3]]);
		/**/
		return leftx;
	}

	/* Set Initial Num of lines -------------------------------------- */
	wn[0]=3;
	/* Save Right part Save for X-Stokes ---------------------*/
	wi[0]=rightx;

	/* Add Coefficients for left parts of X-Stokes ----------------*/
	/* Staggered Nodes num */
	/*                      [0]                [2]   */
	/*                  Ro0,Sxy0,Nu0                 */
	/*                                               */
	/*         Pr1,Sxx1    <Vx0>  Pr3,Sxx3           */
	/*                                               */
	/*                      [1]                [3]   */
	/*                  RO1,Sxy1,Nu1                 */
	/*                                               */
	/*  0(P) 1(Vx)  2(Vy)  */
	/* -dP/dX */
	wn[1]=v[1]*3+0;
	wi[1]=+1.0/xkf;
	wn[2]=v[3]*3+0;
	wi[2]=-1.0/xkf;
	/* DVx/Dt*RO */
	wn[3]=v[0]*3+1;
	wi[3]=-inertyn*mrx[v[0]]/timestepe;
	/* dSIGxx/dX */
	sxxcalc(m1,m2+1,-1.0/xkf);
	sxxcalc(m1+1,m2+1,1.0/xkf);
	/* dSIGxy/dY */
	sxycalc(m1,m2,-1.0/ykf);
	sxycalc(m1,m2+1,1.0/ykf);
	// Next 11 lines added for drunken sailor instability
	/* -gx/dt/2*(vx*dRHO/dx+vy*dRHO/dy) */
	wn[wn[0]+1]=v[0]*3+1;
	wi[wn[0]+1]=-timestep/2.0*GXKOEF*(ro[v[2]]+ro[v[3]]-ro[v[0]-ynumy]-ro[v[1]-ynumy])/(gx[m1+1]-gx[m1-1])/2.0;
	wn[wn[0]+2]=v[0]*3+2;
	wi[wn[0]+2]=-timestep/2.0*GXKOEF*(gx[m1]-gx[m1-1])/(gx[m1+1]-gx[m1-1])/2.0*(ro[v[1]]-ro[v[0]])/(gy[m2+1]-gy[m2]);
	wn[wn[0]+3]=v[1]*3+2;
	wi[wn[0]+3]=-timestep/2.0*GXKOEF*(gx[m1]-gx[m1-1])/(gx[m1+1]-gx[m1-1])/2.0*(ro[v[1]]-ro[v[0]])/(gy[m2+1]-gy[m2]);
	wn[wn[0]+4]=(v[0]-ynumy)*3+2;
	wi[wn[0]+4]=-timestep/2.0*GXKOEF*(gx[m1+1]-gx[m1])/(gx[m1+1]-gx[m1-1])/2.0*(ro[v[1]]-ro[v[0]])/(gy[m2+1]-gy[m2]);
	wn[wn[0]+5]=(v[1]-ynumy)*3+2;
	wi[wn[0]+5]=-timestep/2.0*GXKOEF*(gx[m1+1]-gx[m1])/(gx[m1+1]-gx[m1-1])/2.0*(ro[v[1]]-ro[v[0]])/(gy[m2+1]-gy[m2]);
	wn[0]+=5;

	return 0;
}
/* Left+Right Side or Err for X-Stokes Equation */




/* LEFT+Right Side or Err for Y-Stokes Equation */
/* Stoks equation initial form */
/* dSIGyy/dY + dSIGxy/dX - dP/dY = RO*DVy/Dt - RO*Gy */
double ystokserr(long int m1, long int m2, int ynerr)
	/* m1,m2 - node X,Y number */
	/* ynerr - Err Calc Y(1)/N(0) */
{
	/* Counters */
	int n1;
	long int v[4];
	/* Err Buf */
	double lefty,righty,nueff;
	/* Distances */
	double xkf=gx[m1+1]-gx[m1],ykf=(gy[m2+1]-gy[m2-1])/2.0;

	/* Staggered Nodes num */
	/*                                               */
	/*                Pr2,Syy2                       */
	/*                                               */
	/*         [0]                  [2]              */
	/*    Ro0,Sxy0,Nu0   <Vy0>   Ro2,Sxy2,Nu2        */
	/*                                               */
	/*                Pr3,Syy3                       */
	/*                                               */
	/*         [1]                  [3]              */
	/*                                               */
	v[0]=m1*ynumy+m2;v[1]=v[0]+1;
	v[2]=v[0]+ynumy;v[3]=v[2]+1;

	/* RIGHT part of Y-Stokes */
	/* dSIGyy/dY + dSIGxy/dX - dP/dY = RO*DVy/Dt - RO*Gy */
	righty  = (-GYKOEF-inertyn*mvy[v[0]]/timestepe)*mry[v[0]];

	/* Return val for LSQ err ----------------------------*/
	if(ynerr==1)
	{
		/* LEFT part of Y-Stokes */
		/* dSIGyy/dY + dSIGxy/dX - dP/dY = - RO*Gy */
		/**/
		/* dSIGyy/dY */
		lefty =(-sxx[v[3]]+sxx[v[2]])/ykf;
		/* dSIGxy/dX */
		lefty+=(sxy[v[2]]-sxy[v[0]])/xkf;
		/* -dP/dY */
		lefty-=(pr[v[3]]-pr[v[2]])/ykf;
		/**/
		/* DVy/Dt*RO */
		lefty-=inertyn*vy[v[0]]*mry[v[0]]/timestepe;
		/**/
		/* Next 2 lines added for drunken sailor instability */
		/* -gy/dt/2*(vy*dRHO/dy+vx*dRHO/dx) */
		lefty-=timestep/2.0*GYKOEF*(vy[v[0]]*(ro[v[1]]+ro[v[3]]-ro[v[0]-1]-ro[v[2]-1])/(gy[m2+1]-gy[m2-1])/2.0+((vx[v[0]]+vx[v[2]])*(gy[m2]-gy[m2-1])+(vx[v[0]-1]+vx[v[2]-1])*(gy[m2+1]-gy[m2]))/(gy[m2+1]-gy[m2-1])/2.0*(ro[v[2]]-ro[v[0]])/(gx[m1+1]-gx[m1]));
		/**/
		/* Effective NU calc */
		nueff=MAXV(ABSV(nu[v[0]]),ABSV(nu[v[2]]));
		/**/
		/* Y-STOKS Error */
		lefty=(lefty-righty)/nueff;
		/**/
		/* Min,Max Value of P Save */
		errbuf[10]=MAXV(pr[v[2]],pr[v[3]]);
		errbuf[11]=MINV(pr[v[2]],pr[v[3]]);
		/**/
		return lefty;
	}

	/* Set Initial Num of lines -------------------------------------- */
	wn[0]=3;
	/* Save Right parts Save for Y-Stokes ---------------------*/
	wi[0]=righty;
	/**/
	/* Add Coefficients for left parts of Y-Stokes ----------------*/
	/* Staggered Nodes num */
	/*                                               */
	/*                Pr2,Syy2                       */
	/*                                               */
	/*         [0]                  [2]              */
	/*    Ro0,Sxy0,Nu0   <Vy0>   Ro2,Sxy2,Nu2        */
	/*                                               */
	/*                Pr3,Syy3                       */
	/*                                               */
	/*         [1]                  [3]              */
	/*                                               */
	/*  0(P) 1(Vx)  2(Vy)  */
	/* -dP/dY */
	wn[1]=v[2]*3+0;
	wi[1]=+1.0/ykf;
	wn[2]=v[3]*3+0;
	wi[2]=-1.0/ykf;
	/* DVx/Dt*Ro */
	wn[3]=v[0]*3+2;
	wi[3]=-inertyn*mry[v[0]]/timestepe;
	/* dSIGyy/dY */
	sxxcalc(m1+1,m2,1.0/ykf);
	sxxcalc(m1+1,m2+1,-1.0/ykf);
	/* dSIGxy/dX */
	sxycalc(m1,m2,-1.0/xkf);
	sxycalc(m1+1,m2,1.0/xkf);
	// Next 12 lines added for drunken sailor instability
	/* -gy/dt/2*(vy*dRHO/dy+vx*dRHO/dx) */
	wn[wn[0]+1]=v[0]*3+2;
	wi[wn[0]+1]=-timestep/2.0*GYKOEF*(ro[v[1]]+ro[v[3]]-ro[v[0]-1]-ro[v[2]-1])/(gy[m2+1]-gy[m2-1])/2.0;
	wn[wn[0]+2]=v[0]*3+1;
	wi[wn[0]+2]=-timestep/2.0*GYKOEF*(gy[m2]-gy[m2-1])/(gy[m2+1]-gy[m2-1])/2.0*(ro[v[2]]-ro[v[0]])/(gx[m1+1]-gx[m1]);
	wn[wn[0]+3]=v[2]*3+1;
	wi[wn[0]+3]=-timestep/2.0*GYKOEF*(gy[m2]-gy[m2-1])/(gy[m2+1]-gy[m2-1])/2.0*(ro[v[2]]-ro[v[0]])/(gx[m1+1]-gx[m1]);
	wn[wn[0]+4]=(v[0]-1)*3+1;
	wi[wn[0]+4]=-timestep/2.0*GYKOEF*(gy[m2+1]-gy[m2])/(gy[m2+1]-gy[m2-1])/2.0*(ro[v[2]]-ro[v[0]])/(gx[m1+1]-gx[m1]);
	wn[wn[0]+5]=(v[2]-1)*3+1;
	wi[wn[0]+5]=-timestep/2.0*GYKOEF*(gy[m2+1]-gy[m2])/(gy[m2+1]-gy[m2-1])/2.0*(ro[v[2]]-ro[v[0]])/(gx[m1+1]-gx[m1]);
	wn[0]+=5;

	return 0;
}
/* Left+Right Side or Err for Y-Stokes Equation */




/* Left side or Err for vX Boundary Condition Equation */ 
/* V=CONST+KOEF*Vn */
double xbonderr(long int mcmax, int ynerr)
	/* mcmax - numer of cur Vx in sol[] */
	/* ynerr - Err Calc Y(1)/N(0) */
{
	/* Val Buffer */
	double leftx=0;
	int n1;

	/* Error Calc */
	if (ynerr)
	{
		/* Add Const */
		leftx=x[mcmax]-bondv[bondm[mcmax]][0];
		/* Add Koef */
		for (n1=0;n1<3;n1++)
		{
			if(bondn[bondm[mcmax]][n1]) 
			{
				leftx-=bondv[bondm[mcmax]][n1+1]*x[bondn[bondm[mcmax]][n1]-1];
			}
		}
		/**/
		return leftx;
	}
	/* Add X CONST */
	wn[0]=1;
	wi[0]=bondv[bondm[mcmax]][0];
	
	/* ============================================== */	
        /* Convergent velocity reduction during collision */
	timecc0 = 10.0e+6;
	timecc1 = 100.0e+6;
	/* ============= */	
	vel_cc2 = 0.0000;     /* 0.0 cm/yr */
	vel_cc1 = 0.0000;     /* 0.0 cm/yr */
	vel_cc0 = MAXV(vel_cc0,wi[0]);

        if (collision)
        {
	  for (n1=0;n1<=wn[0];n1++)
              {
	        if ((timesum/3.15576e+7)>timecc1)
	           {	
	            if (wi[n1] > 0)  {wi[n1]= vel_cc1;}
	            if (wi[n1] < 0)  {wi[n1]=-vel_cc2;}
	           }
                if  ((timesum/3.15576e+7)>timecc0 && (timesum/3.15576e+7)<timecc1)
                   {
		    vel_cc = vel_cc0-(vel_cc0-vel_cc1)*((timesum/3.15576e+7)-timecc0)/(timecc1-timecc0);
                    if (wi[n1] > 0) wi[n1]= vel_cc;
		    vel_cc = vel_cc0-(vel_cc0-vel_cc2)*((timesum/3.15576e+7)-timecc0)/(timecc1-timecc0);
		    if (wi[n1] < 0) wi[n1]=-vel_cc;
		   }
	      }
        }
	//if (printmod)  {fprintf(fp_log," Vx = %e \n ",wi[0]); fflush(fp_log);}
	/* === end === */

	wn[1]=mcmax;
	wi[1]=1.0;
	/* Add X PAR1,PAR2,PAR3 */
	for (n1=0;n1<3;n1++)
	{
		if(bondn[bondm[mcmax]][n1]) 
		{
			wn[0]+=1;
			wn[wn[0]]=bondn[bondm[mcmax]][n1]-1;
			wi[wn[0]]=-bondv[bondm[mcmax]][n1+1];
		}
	}
	return 0;
}
/* Left side or Err for vX Boundary Condition Equation */ 





/* Left side or Err for vY Boundary Condition Equation */ 
/* V=CONST+KOEF*Vn */
double ybonderr(long int mcmax, int ynerr)
	/* mcmax - numer of cur Vx in sol[] */
	/* ynerr - Err Calc Y(1)/N(0) */
{
	/* Val Buffer */
	double lefty=0;
	int n1;

	/* Error Calc */
	if (ynerr)
	{
		/* Add Const */
		lefty=x[mcmax]-bondv[bondm[mcmax]][0];
		/* Add Koef */
		for (n1=0;n1<3;n1++)
		{
			if(bondn[bondm[mcmax]][n1]) 
			{
				lefty-=bondv[bondm[mcmax]][n1+1]*x[bondn[bondm[mcmax]][n1]-1];
			}
		}
		/**/
		return lefty;
	}

	/* Add Y CONST */
	wn[0]=1;
	wi[0]=bondv[bondm[mcmax]][0];
	wn[1]=mcmax;
	wi[1]=1.0;
	/* Add Y PAR1,PAR2,PAR3 */
	for (n1=0;n1<3;n1++)
	{
		if(bondn[bondm[mcmax]][n1]) 
		{
			wn[0]+=1;
			wn[wn[0]]=bondn[bondm[mcmax]][n1]-1;
			wi[wn[0]]=-bondv[bondm[mcmax]][n1+1];
		}
	}

	return 0;
}
/* Left side or Err for vY Boundary Condition Equation */ 




/* Left side or Err for P Boundary Condition Equation */ 
/* P=CONST+KOEF*Pn */
double pbonderr(long int mcmax, int ynerr)
	/* mcmax - numer of cur P  in sol[] */
	/* ynerr - Err Calc Y(1)/N(0) */
{
	/* Val Buffer */
	double leftp;
	int n1;

	/* Error Calc */
	if (ynerr)
	{
		/* Add Const */
		leftp=x[mcmax]-bondv[bondm[mcmax]][0];
		/* Add Koef */
		for (n1=0;n1<3;n1++)
		{
			if(bondn[bondm[mcmax]][n1]) 
			{
				leftp-=bondv[bondm[mcmax]][n1+1]*x[bondn[bondm[mcmax]][n1]-1];
			}
		}
		return leftp*leftp;
	}

	/* Add P CONST */
	wn[0]=1;
	wi[0]=bondv[bondm[mcmax]][0];
	wn[1]=mcmax;
	wi[1]=1.0;
	/* Add P PAR1,PAR2,PAR3 */
	for (n1=0;n1<3;n1++)
	{
		if(bondn[bondm[mcmax]][n1]) 
		{
			wn[0]+=1;
			wn[wn[0]]=bondn[bondm[mcmax]][n1]-1;
			wi[wn[0]]=-bondv[bondm[mcmax]][n1+1];
		}
	}
	return 0;
}
/* Left side or Err for P Boundary Condition Equation */ 




/* Weight of FD calculation for after Fornberg (1996) */
void fdweight(int n, int m, double xi)
	/* n - maximal index 0-n */
	/* m - required derivative order 0-m */
	/* xi - derivation point coordinate */
{
	/* Counters */
	int i,j,k,mn;
	double c1,c2,c3,c4,c5,kk;

	c1=1.0;
	c4=xn[0]-xi;
	for(k=0;k<=m;k++)
	{
		for(j=0;j<=n;j++)
		{
			cn[j][k]=0;
		}
	}

	cn[0][0]=1.0;
	for(i=1;i<=n;i++)
	{
		mn=i;if(mn>m) mn=m;
		c2=1.0;
		c5=c4;
		c4=xn[i]-xi;
		for(j=0;j<i;j++)
		{
			c3=xn[i]-xn[j];
			c2*=c3;
			for(k=mn;k>0;k--)
			{
				kk=(double)(k);
				cn[i][k]=c1*(kk*cn[i-1][k-1]-c5*cn[i-1][k])/c2;
			}
			cn[i][0]=-c1*c5*cn[i-1][0]/c2;
			for(k=mn;k>0;k--)
			{
				kk=(double)(k);
				cn[j][k]=(c4*cn[j][k]-kk*cn[j][k-1])/c3;
			}
			cn[j][0]=c4*cn[j][0]/c3;
		}
		c1=c2;
	}
}
/* Weight of FD calculation after Fornberg (1996) */



/* Left side or Value for Sxx  Equation */ 
/* Sxx=2Nu*Exx*Z+Sxx0*(1-Z), Exx=1/2(dVx/dX-dVy/dY) */
double sxxcalc(long int m1, long int m2, double ynval)
	/* m1,m2 - node X,Y number */
	/* ynval - Val Sxx Calc Y(0)/N(koefficient) */
{
	/* Exx horizontal position */
	double xi=(gx[m1-1]+gx[m1])/2.0,leftsxx=0,nueff,sxxeeff,xelvis=1.0,ggeff=0;
	long int v[4];
	int n1,n;

	/* Staggered Nodes num */
	/*   [0]      Vy0       [2] */
	/*   Nu0                Nu2 */
	/*                          */
	/*   Vx0    Sxx3,Exx3   Vx2 */
	/*                          */
	/*   [1]                [3] */
	/*   Nu1      Vy1       Nu3 */
	/*                          */
	v[0]=(m1-1)*ynumy+(m2-1);v[1]=v[0]+1;
	v[2]=v[0]+ynumy;v[3]=v[2]+1;

	/* Effective viscosity calc */
	nueff=nd[v[3]];

	/* Effective elastic stress calc */
	// Calculate Sxx with visco-elasticity factor following eq. 13 in Gerya & Yuen (2007)
	if(stoksmod)
	{
		/* Effective shear modulus calc */
		ggeff=gd[v[3]];
		/* Effective viscoelastic factor calc */
		xelvis=ggeff*timestepe/(ggeff*timestepe+nueff);
		/* Effective old elastic stress calc, save */
		sxxeeff=sxxe[v[3]];
		/* Return viscoelasticity factor */
		if(stoksmod<0)
		{
			eps[0]=nueff;
			eps[1]=ggeff;
			return xelvis;
		}

		/* Nu recalc, Check */
		nueff*=xelvis; if(nueff<nubeg) nueff=nubeg;
	}
	/* Nu Save */
	eps[2]=nueff;

	/* Staggered Nodes num */
	/*   [0]      Vy0       [2] */
	/*   Nu0                Nu2 */
	/*                          */
	/*   Vx0    Sxx3,Exx3   Vx2 */
	/*                          */
	/*   [1]                [3] */
	/*   Nu1      Vy1       Nu3 */
	/*                          */

	/* Return Sxx,Exx val ----------------------------*/
	// The times xelvis in the first mainly viscous part of the equation is included in nueff
	if(ynval==0)
	{
		/* Exx=1/2(dVx/dX-dVy/dY) */
		leftsxx=0.5*((vx[v[2]]-vx[v[0]])/(gx[m1]-gx[m1-1])-(vy[v[1]]-vy[v[0]])/(gy[m2]-gy[m2-1]));

		/* Save Exx */
		eps[0]=leftsxx;

		/* Calc Sxx=2Nu*Exx*X+Sxx0*(1-X) */
		leftsxx=2.0*nueff*leftsxx;
		/* Effective elastic stress calc */
		if(stoksmod)
		{
			leftsxx+=sxxeeff*(1.0-xelvis);
		}

		return leftsxx;
	}

	/* Add Coefficients for left parts of Sxx ----------------*/
	/*  0(P) 1(Vx)  2(Vy)  */
	/* Sxx=2Nu*Exx*Z+Sxx0*(1-Z), Exx=1/2(dVx/dX-dVy/dY) */
	/* Add Vx with koefficients */
	wn[wn[0]+1]=v[0]*3+1;
	wi[wn[0]+1]=-ynval*nueff/(gx[m1]-gx[m1-1]);
	wn[wn[0]+2]=v[2]*3+1;
	wi[wn[0]+2]=+ynval*nueff/(gx[m1]-gx[m1-1]);
	/* Add Vy with koefficients */
	wn[wn[0]+3]=v[0]*3+2;
	wi[wn[0]+3]=+ynval*nueff/(gy[m2]-gy[m2-1]);
	wn[wn[0]+4]=v[1]*3+2;
	wi[wn[0]+4]=-ynval*nueff/(gy[m2]-gy[m2-1]);

	/* Add elastic stress to the right part */
	if(stoksmod)
		{wi[0]-=(1.0-xelvis)*ynval*sxxeeff;}

	/* Add total number of lines */
	wn[0]+=4;

	return 0;
}
/* Left side or Value for Sxx  Equation */ 




/* Left side or Value for Sxy  Equation */ 
/* Sxy=2Nv*Exy, Exy=1/2(dVx/dY+dVy/dX) */
double sxycalc(long int m1, long int m2, double ynval)
	/* m1,m2 - node X,Y number */
	/* ynval - Val Syy Calc Y(0)/N(koefficient) */
{
	/* Exy position */
	double xi,leftsxy=0,leftsxy1=0,leftsxy2=0,nueff,sxyeeff,xelvis=1.0,ggeff=0;
	long int v[4];
	int n1;

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
	v[0]=(m1-1)*ynumy+(m2-1);v[1]=v[0]+1;
	v[2]=v[0]+ynumy;v[3]=v[2]+1;

	/* Effective viscosity calc */
	nueff=nu[v[3]];

	/* Effective elastic stress calc */
	if(stoksmod)
	{
		/* Effective shear modulus calc */
		ggeff=gg[v[3]];
		/* Effective viscoelastic factor calc */
		xelvis=ggeff*timestepe/(ggeff*timestepe+nueff);
		/* Effective old elastic stress calc, save */
		sxyeeff=sxye[v[3]];
		/* Return viscoelasticity factor */
		if(stoksmod<0)
		{
			eps[0]=nueff;
			eps[1]=ggeff;
			return xelvis;
		}

		/* Nu recalc, Check */
		nueff*=xelvis; if(nueff<nubeg) nueff=nubeg;
	}
	/* Nu Save */
	eps[2]=nueff;

	/* Staggered Nodes num */
	/*  [0]                [2]            */
	/*                                    */
	/*                     Vx2            */
	/*                                    */
	/*  [1]     Vy1        [3]       Vy3  */
	/*                   Exy3,Nu3         */
	/*                                    */
	/*                     Vx3            */
	/* Return Sxy,Exy val ----------------------------*/
	if(ynval==0)
	{
		/* Exy=1/2(dVx/dY+dVy/dX)=0 */
		leftsxy1=(vx[v[3]]-vx[v[2]])/(gy[m2+1]-gy[m2-1]);
		leftsxy2=(vy[v[3]]-vy[v[1]])/(gx[m1+1]-gx[m1-1]);

		/* Save Exy */
		eps[0]=leftsxy=leftsxy1+leftsxy2;
		/* Save Esp (rotation rate) */
		eps[1]=leftsxy1-leftsxy2;

		/* Calc Sxy=2Nu*Exy-Sxyp */
		leftsxy=2.0*nueff*leftsxy;
		/* Effective elastic stress calc */
		if(stoksmod)
		{
			leftsxy+=sxyeeff*(1.0-xelvis);
		}
		return leftsxy;
	}

	/* Add Coefficients for left parts of Sxy ----------------*/
	/*  0(P) 1(Vx)  2(Vy)  */
	/* Sxy=2Nu*Exy, Exy=1/2(dVx/dY+dVy/dX) */
	/* Add Vx with koefficients */
	wn[wn[0]+1]=v[2]*3+1;
	wi[wn[0]+1]=-ynval*2.0*nueff/(gy[m2+1]-gy[m2-1]);
	wn[wn[0]+2]=v[3]*3+1;
	wi[wn[0]+2]=+ynval*2.0*nueff/(gy[m2+1]-gy[m2-1]);
	/* Add Vy with koefficients */
	wn[wn[0]+3]=v[1]*3+2;
	wi[wn[0]+3]=-ynval*2.0*nueff/(gx[m1+1]-gx[m1-1]);
	wn[wn[0]+4]=v[3]*3+2;
	wi[wn[0]+4]=+ynval*2.0*nueff/(gx[m1+1]-gx[m1-1]);

	/* Add elastic stress to the right part */
	if(stoksmod)
	{
		wi[0]-=(1.0-xelvis)*ynval*sxyeeff;
	}

	/* Add total Num of lines */
	wn[0]+=4;

	return 0;
}
/* Left side or Value for Sxy  Equation */ 




/* Left side or Err for Compressible Continuity Equation  */
/* div(V) = -D(ln(RO))/dt, div(V)=dVx/dX+dVy/dY */
double conterr(long int m1, long int m2, int ynerr)
	/* m1,m2 - node X,Y number */
	/* ynerr - Err Calc Y(1)/N(0) */
{
	/* Counter */
	long int v[4];
	/* Val Buffer */
	double leftc=0,rightc=0,dlnrodp;
	int n1;

	/* Staggered Nodes num */
	/*   [0]       Vy0      [2] */
	/*                          */
	/*   Vx0        <P3>    Vx2 */
	/*            Exx3,Eyy3     */
	/*                          */
	/*   [1]       Vy1      [3] */
	v[0]=(m1-1)*ynumy+(m2-1);v[1]=v[0]+1;
	v[2]=v[0]+ynumy;v[3]=v[2]+1;
	/**/
	/**/
	/**/
	/* -D(ln(RO))/dt add to the right part */
	rightc=0;
	/* D(ln(RO))/dp*dP/dt add to the right part */
	dlnrodp=drp[v[3]];
	if(timestepe>0)
	{
		rightc=dro[v[3]];
		rightc=(dlnrodp*sppe[v[3]]-rightc)/timestepe;
	}

	/* Staggered Nodes num */
	/*   [0]       Vy0      [2] */
	/*                          */
	/*   Vx0        <P3>    Vx2 */
	/*            Exx3,Eyy3     */
	/*                          */
	/*   [1]       Vy1      [3] */
	/* Return dVx/dX+dVy/dY err ----------------------------*/
	if(ynerr==1)
	{
		/* div(V)=dVx/dX+dVy/dY */
		leftc=(vx[v[2]]-vx[v[0]])/(gx[m1]-gx[m1-1])+(vy[v[1]]-vy[v[0]])/(gy[m2]-gy[m2-1]);
		/* D(ln(RO))/dp*dP/dt add to the left part */
		if(timestepe>0 && dlnrodp)
		{
			leftc+=dlnrodp*pr[v[3]]/timestepe;
		}
		return leftc-rightc;
	}

	/* Add continuity equation */
	/* Set Initial Num of lines -------------------------------------- */
	wn[0]=0;

	/* Save Right part for Contin ---------------------*/
	wi[0]=rightc;

	/* Add Coefficients for left parts of div(V) ----------------*/
	/*  0(P) 1(Vx)  2(Vy)  */
	/* div(V)=dVx/dX+dVy/dY */
	/* Add Vx with koefficients */
	wn[1]=v[0]*3+1;
	wi[1]=-1.0/(gx[m1]-gx[m1-1]);
	wn[2]=v[2]*3+1;
	wi[2]=+1.0/(gx[m1]-gx[m1-1]);
	/* Add Vy with koefficients */
	wn[3]=v[0]*3+2;
	wi[3]=-1.0/(gy[m2]-gy[m2-1]);
	wn[4]=v[1]*3+2;
	wi[4]=+1.0/(gy[m2]-gy[m2-1]);
	/* Add total Num of lines */
	wn[0]=4;

	/* D(ln(RO))/dp*dP/dt add to the left part */
	if(timestepe>0 && dlnrodp)
	{
		wn[5]=v[3]*3;
		wi[5]=dlnrodp/timestepe;
		wn[0]=5;
	}

	/* Check Boundary conditions around Cell */
	leftc=1.0;
	if (!bondm[v[0]*3+1]) leftc=0;
	if (!bondm[v[0]*3+2]) leftc=0;
	if (!bondm[v[1]*3+2]) leftc=0;
	if (!bondm[v[2]*3+1]) leftc=0;

	return leftc;
}
/* Left side or Err for Continuity Equation */

