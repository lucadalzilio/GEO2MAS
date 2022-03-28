/* OMP Calculation of new T after time step */
void titerateomp(int m0)
	/* m0 - Circle number */
{
	/* Counters */
	long int m1,m2,m3,m4,mcmax,mcmax1;
	int n1;
	/**/
	/* Val buffer */
	double ival,maxtkrate,tktimestep=0,tktimesum,tktimesum0,celdx,celdy,swt1,maxtk,mintk;
	/* Err koef */
	double heatsum=0,heatnum=0,bondsum=0,bondnum=0,mrocp,mkt,dtk,dtk1;
	// Parallelization
	double start,TK,TK2;
	long int m10,m20;
	int nt;

#pragma omp parallel 
	{
		nt=omp_get_num_threads();
	}
	
	/* Right par of equation T step limit definition */
	/* Save Old temperatures */
	/* Save TK */
	for (m1=0;m1<nodenum;m1++)
	{
		tk0[m1]=tk[m1];
	}
	/**/
	maxtkrate=0;
	for (m1=1;m1<xnumx-1;m1++)
		for (m2=1;m2<ynumy-1;m2++)
	{
		/* Cur Line Num in sol0[] */
		m3=nodenum3+m1*ynumy+m2;
		if(bondm[m3]==0)
		{
			/* Max dT/dt definition */
			ival=heaterr(m1,m2,0);
			ival=ABSV(ival);
			maxtkrate=MAXV(maxtkrate,ival);
		}
	}

	/* Time step calc */
	if (maxtkrate)
	{
		tktimestep=maxtkstep/maxtkrate;
		if (printmod) fprintf(fp_log,"\n !!! MAX VALID TIME STEP FOR TEMPERATURE %e YEAR !!!\n",tktimestep/3.15576e+7);
		fflush(fp_log);
	}

	/* Check timestep */
	if (timestep==0 || tktimestep==0) return;

	/* Thermal Solution Cycle ------------------------------------------ */
	/* Save TK */
	for (m1=0;m1<nodenum;m1++)
	{
		tk2[m1]=tk[m1];
	}

	tktimesum=0;
	tktimesum0=timestep;
	do
	{
		/* Timestep for temperature selection */
		timestep=MINV(tktimesum0,tktimestep);
		if (timestep>(tktimesum0-tktimesum)) timestep=tktimesum0-tktimesum;

		/* Clear solution */
		pos0cur=0;
		/* Save TK */
		for (m1=0;m1<nodenum;m1++)
		{
			tk0[m1]=tk[m1];
		}

		/* Add Matrix by TK Equations */
		if (printmod) fprintf(fp_log,"Adding Matrix By TK ...");
		fflush(fp_log);
		/* TK - Node Cycle */
		for (m1=0;m1<xnumx;m1++)
			for (m2=0;m2<ynumy;m2++)
		{
			/* Cur Line Num in sol0[] */
			mcmax=m1*ynumy+m2;
			mcmax1=mcmax+nodenum3;
			if(!bondm[mcmax1])
			{
				/* Right part calculation */
				heaterr(m1,m2,0);
				/* Thermal conductivity Equation koef calc add */
				heaterr(m1,m2,1);

				/* Add matrix */
				gausmat4(3,mcmax,0);
			}
			else	
			{
				/* Add vX Simple Boundary */
				tbonderr(mcmax,0);

				/* Rescale coefficients */
				for (n1=0;n1<=wn[0];n1++) wi[n1]/=timestep;
	
				/* Add matrix */
				gausmat4(3,mcmax,0);
			}
		}
		if (printmod) fprintf(fp_log,"OK!\n");
		fflush(fp_log);
		/* End  Add Matrix By TK Equations */

		/* Solve Matrix */
		start=omp_get_wtime();
		gausmat4(0,nodenum,0);
		if (printmod == 10000) fprintf(fp_log," Time taken for solving temperature equations = %e s ",omp_get_wtime()-start);
		/* End Solve Matrix */

		/* TK reload */
		maxtk=-1e+30;
		mintk=1e+30;
		for (m1=0;m1<nodenum;m1++)
		{
			tk[m1]=x[m1];
			maxtk=MAXV(maxtk,tk[m1]);
			mintk=MINV(mintk,tk[m1]);
		}

		/* Error check cycle */
		heatsum=0;
		heatnum=0;
		for (m1=0;m1<xnumx;m1++)
			for (m2=0;m2<ynumy;m2++)
		{
			/* Cur Line Num in sol0[] */
			mcmax=m1*ynumy+m2;
			mcmax1=mcmax+nodenum3;
			if(!bondm[mcmax1])
			{
				/* Heat equat Err */
				ival=heaterr(m1,m2,2);
				heatsum+=ival*ival;
				heatnum+=1.0;
				/* Print Results */
				ival/=(kfxx+kfyy)*0.5;
				if (printmod && ABSV(ival)>HEATMIN)
				{
					fprintf(fp_log,"\n Large T K  err at X=%ld Y=%ld:   Err=%e",m1,m2,ival);fflush(fp_log);
				}
			}
			else
			{
				/* Add vX-vY-Boundary */
				ival=tbonderr(mcmax,1);
				bondsum+=ival*ival;
				bondnum+=1.0;
			}
		}
		heatsum=pow(heatsum/heatnum,0.5)/(kfxx+kfyy)*2.0;
		bondsum=pow(bondsum/bondnum,0.5);

		/* Add thermal step */
		tktimesum+=timestep;

		/* Print Results */
		if (printmod)
		{
			fprintf(fp_log,"\n KRUG %2d \n",m0+1);
			fprintf(fp_log,"THERMAL STEP = %e YEARS    THERMAL TIME = %e YEARS \n",timestep/3.15576e+7,tktimesum/3.15576e+7);
			fprintf(fp_log,"TEMPERATURE MIN = %e MAX = %e \n",mintk,maxtk);
			fprintf(fp_log,"TERMAL EQUATION: num = %e err = %e \n",heatnum,heatsum);
			fprintf(fp_log,"BOUND T: num = %e err = %e \n",bondnum,bondsum);
			fflush(fp_log);
		}
	}
	while(tktimesum<tktimesum0);
	/* Restore timestep */
	timestep=tktimesum0;
	/* End Thermal Solution Cycle ------------------------------------------ */

	/* --- Recalc temperature for markers --- */
	
	/* Clear nodes wt, save increment */
	for (m1=0;m1<nodenum;m1++)
		{tol0[m1]=tol1[m1]=0;}
	
	start=omp_get_wtime();
		
#pragma omp parallel shared(marknum,markx,marky,markt,markk,xsize,ysize,xnumx,ynumy,gx,gy,tk,tk2,tdeep,heatdif,kt,cp,ro,timestep) \
	private(TK,TK2,dtk,dtk1,celdx,celdy,swt1,mkt,mrocp,ival)
	{	
		// Initialize temporarily interpolation arrays (capitilized) to zero (inside pragma, so is private)
		double* Tol0 	= (double*) calloc(nodenum,sizeof(double));
		double* Tol1 	= (double*) calloc(nodenum,sizeof(double));
	
		/* Reset T for fluid markers */
#pragma omp for schedule(runtime)
		for (int m3=0;m3<marknum;m3++) 
		{
			/* Check markers out of grid */
			if (markx[m3]>0 && marky[m3]>0 && markx[m3]<xsize && marky[m3]<ysize && markt[m3]>=50 && markt[m3]<100)
			{
				int m10=m1serch(markx[m3]);
				int m20=m2serch(marky[m3]);
			
				/* Interpolate temperature */
				allintertomp(markx[m3],marky[m3],m10,m20,&TK,&TK2);
				markk[m3]=TK;
			}
		}

		/* Recalc marker T + Diffusion, Add  nodes wt */
	
		/* Define step */
#pragma omp for schedule(runtime)
		for (int m3=0;m3<=marknum;m3++) 
		{
			// For normal and molten rock types
			/* Check markers out of grid */
			if (markx[m3]>0 && marky[m3]>0 && markx[m3]<xsize && marky[m3]<ysize && markt[m3]<50)
			{
				int m10=m1serch(markx[m3]);
				int m20=m2serch(marky[m3]);
			
				/* Interpolate temperature */
				allintertomp(markx[m3],marky[m3],m10,m20,&TK,&TK2);

				/* Reset marker temperature for newly coming markers */
				if(markk[m3]<=0) 
				{
					markk[m3]=TK2;
					if(markk[m3]<tdeep) markk[m3]=tdeep;
				}

				/* Numerical diffusion add to markers -------------*/
				if(heatdif)
				{
					/* Calc difference in marker temperature */
					dtk=TK2-markk[m3];

					/* Interpolation of k, ro, Cp from nodes to marker */
					/* Marker weight calculation using dimension of current Cell */
					celdx=gx[m10+1]-gx[m10];
					celdy=gy[m20+1]-gy[m20];
					swt1=1.0/celdx/celdy;

					celdx=(markx[m3]-gx[m10])/(gx[m10+1]-gx[m10]);
					celdy=(marky[m3]-gy[m20])/(gy[m20+1]-gy[m20]);
					if (celdx<0 || celdy<0 || celdx>1.0 ||celdy>1.0 || swt1<0) {fprintf(fp_log,"ERROR..: num=%ld x=%e y=%e celdx=%e celdy=%e swt1=%e",m3,markx[m3],marky[m3],celdx,celdy,swt1);fflush(fp_log);getchar();}
				
					/* Interpolate k, ro, Cp using interpolation koefficients */
					mkt=mrocp=0;
				
					/* Reload horizontal coefficients to cn[] */
					ival = 0;
					for (int m1=m10;m1<=m10+1;m1++)
					{
						for (int m2=m20;m2<=m20+1;m2++)
						{
							/* Current node num, wt */
							int m4=m1*ynumy+m2;

							if (m1==m10 && m2==m20)
							{
								ival=(1.0-celdx)*(1.0-celdy);
							}
							else if (m1==m10+1 && m2==m20)
							{
								ival=(celdx)*(1.0-celdy);
							}
							else if (m1==m10 && m2==m20+1)
							{
								ival=(1.0-celdx)*(celdy);
							}
							else if (m1==m10+1 && m2==m20+1)
							{
								ival=(celdx)*(celdy);
							}
						
							mkt+=kt[m4]*ival;
							mrocp+=cp[m4]*ro[m4]*ival;
						}
					}
				
					/* Calc, check diffusion changes in marker temperature: dT/dt=k/ro/cp*(d2T/dX2+d2T/dY2) */
					/* dT=k/ro/cp*(Tgrid-2Tmark+Tgrid)*(1/(xstep)^2+1/(ystep)^2)*dt */
					celdx=gx[m10+1]-gx[m10];
					celdy=gy[m20+1]-gy[m20];
				
					dtk1=-heatdif*mkt/mrocp*2.0*(1.0/celdx/celdx+1.0/celdy/celdy)*timestep;
					if(dtk1<-150.0) dtk1=-150.0;
					dtk1=dtk*(1.0-exp(dtk1));

					/* Diffuse Marker Temperature */
					markk[m3]+=dtk1;
				
					/* Wt for nodes calc, add */
					// Redo, since overwrote celdx,celdy, though swt1 from start this loop & if is still valid
					celdx=(markx[m3]-gx[m10])/(gx[m10+1]-gx[m10]);
					celdy=(marky[m3]-gy[m20])/(gy[m20+1]-gy[m20]);
					if (celdx<0 || celdx>1 || celdy <0 || celdy>1 || swt1<0 ) {fprintf(fp_log,"ERROR..: %ld %e %e %e \n",m3,celdx,celdy,swt1); fflush(fp_log); getchar();}
				
					for (int m1=m10;m1<=m10;m1++)
					{
						for (int m2=m20;m2<=m20;m2++)
						{
							/* Cur Node Num, wt */
							int m4=m1*ynumy+m2;
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
							Tol0[m4]+=ival;
							Tol1[m4]+=dtk1*ival;
						}
					}
				}
				/* End Numerical diffusion add to markers -------------*/

				/* Change marker temperature after solution if no diffusion */
				markk[m3]+=(TK-TK2);
			}
		}
	
		// Add interpolation arrays from different processors and free their memory
#pragma omp critical (sumtolarrays)
		{
			for (int m5=0;m5<nodenum;m5++)
			{
				tol0[m5]+=Tol0[m5];
				tol1[m5]+=Tol1[m5];	   
			}	
		}
		
		// Free dynamically allocated interpolation arrays	
		free(Tol0);
		free(Tol1); 
	}
	// End OMP loop add diffusion T

	/* Numerical antidiffusion add to markers -------------*/
	if(heatdif)
	{
		/* Recalc changes in nodes T */
		for (m1=0;m1<nodenum;m1++)
		{
			if(tol0[m1]) 
			{
				/* Averaged Changes in temperature due to smoothing */
				tol1[m1]/=tol0[m1];
			}
		}

		/* Recalculate marker T + Antidiffusion for normal and molten rock markers */
#pragma omp parallel for shared(marknum,markx,marky,markt,markk,tol1,xsize,ysize,xnumx,ynumy,gx,gy,tk,tk2) \
		private(m1,m2,m3,m4,m10,m20,celdx,celdy,ival) \
			schedule(runtime)
				for (m3=0;m3<=marknum;m3++) 
		{
			/* Check markers out of grid */
			if (markx[m3]>0 && marky[m3]>0 && markx[m3]<xsize && marky[m3]<ysize && markt[m3]<50)
			{
				m10=m1serch(markx[m3]);
				m20=m2serch(marky[m3]);
				
				/* Wt for nodes calc, add */
				celdx=(markx[m3]-gx[m10])/(gx[m10+1]-gx[m10]);
				celdy=(marky[m3]-gy[m20])/(gy[m20+1]-gy[m20]);
				if (celdx<0 || celdx>1 || celdy <0 || celdy>1) {fprintf(fp_log,"ERROR: %ld %e %e %e \n",m3,celdx,celdy,swt1); fflush(fp_log); getchar();}
				
				for (m1=m10;m1<=m10;m1++)
				{
					for (m2=m20;m2<=m20;m2++)
					{
						/* Cur Node Num, wt */
						m4=m1*ynumy+m2;

						if (m1==m10 && m2==m20)
						{
							ival=(1.0-celdx)*(1.0-celdy);
						}
						else if (m1==m10+1 && m2==m20)
						{
							ival=(celdx)*(1.0-celdy);
						}
						else if (m1==m10 && m2==m20+1)
						{
							ival=(1.0-celdx)*(celdy);
						}
						else if (m1==m10+1 && m2==m20+1)
						{
							ival=(celdx)*(celdy);
						}
						
						/* Antidiffuse Marker Temperature */
						if(tol1[m4])
						{
							markk[m3]-=(tol1[m4]*ival);
						}
					}
				}
			}
		}
		// End OMP anti diffusion loop
	}
	/* End Numerical antidiffusion add to markers -------------*/

	if (printmod==10000) fprintf(fp_log,"Time for temperature (anti-)diffusion = %e s\n",omp_get_wtime()-start);
	
	/* End T interpolation ------------------------ */
}
/* End OMP calculation of new T after time step */



/* Left side, Right side or Err of Lagrangian Thermal Conductivity equation calc */
/* Explicit form  */
/* dT/dt = 1/Cp/ro*[dQx/dX + dQy/dY + dH/dt] */
/* Implicit form */
/* dT/dt - 1/Cp/ro*[dQx/dX + dQy/dY] =  1/Cp/ro*(dH/dt] */
double heaterr(long int m1, long int m2, int ynerr)
	/* m1,m2 - node X,Y number */
	/* ynerr - Calc mode: 0-Right, 1-Left, 2-Err */
{
	/* Counters */
	int n1;
	long int v[9];
	/* Buffer */
	double mpb,ival=0,ival1=0,ival2=0,epsval,dpdx,dpdy,e,n;
	/* Distances */
	double xkf=(gx[m1+1]-gx[m1-1])/2.0,ykf=(gy[m2+1]-gy[m2-1])/2.0;

	/* T-Nodes num */
	/* 0  3  6 */
	/* 1  4  7 */
	/* 2  5  8 */
	v[0]=(m1-1)*ynumy+(m2-1);v[1]=v[0]+1;v[2]=v[1]+1;
	v[3]=v[0]+ynumy;v[4]=v[3]+1;v[5]=v[4]+1;
	v[6]=v[3]+ynumy;v[7]=v[6]+1;v[8]=v[7]+1;

	/* Mode of Calculation */
	switch(ynerr)
	{
		/* Right Side Calculation ---------------------- */
		case 0:
		
		/* Vx, Vy, EPS,SIG, in current node calc   */
		allinters(gx[m1],gy[m2]);

		/* Lagrangian Thermal conductivity Equation in Simple Direct form use  */
		/* dT/dt = K/Cp/ro*(d2T/dX2+d2T/dY2)+1/Cp/ro*(dT/dX*dK/dX+dT/dY*dK/dY) + dH/dt/Cp/ro */
		/**/
		/* Heat Sources dH/dt/Cp/ro */
		/* Radioactive heat dH/dt, Wt/m3 */
		ival+=ht[v[4]]/cp[v[4]]/ro[v[4]];
		/**/
		/* Correct adiabatic term alphV*T*(vX*dP/dX+vY*dP/dY)/Cp/ro */
		if(adiabyn==2)
		{
			e=(gx[m1]-gx[m1-1])/(gx[m1+1]-gx[m1-1]);
			n=(gy[m2]-gy[m2-1])/(gy[m2+1]-gy[m2-1]);
			dpdx=2.0*((1.0-n)*(pr[v[7]]-pr[v[4]])+n*(pr[v[8]]-pr[v[5]]))/(gx[m1+1]-gx[m1-1]);
			dpdy=2.0*((1.0-e)*(pr[v[5]]-pr[v[4]])+e*(pr[v[8]]-pr[v[7]]))/(gy[m2+1]-gy[m2-1]);
			ival+=et[v[4]]*tk0[v[4]]/ro[v[4]]/cp[v[4]]*(eps[11]*dpdx+eps[12]*dpdy);
			/* Adiabate computing */
			if(1==0 && timesum<(3.15576e+7*1e+3))
			{
				mpb=(1.0-e)*(1.0-n)*pr[v[4]]+(1.0-e)*n*pr[v[5]]+e*(1.0-n)*pr[v[7]]+e*n*pr[v[8]];
				ival+=et[v[4]]*tk0[v[4]]/ro[v[4]]/cp[v[4]]*mpb/(3.15576e+7*1e+3);
			}
			/*
			ival1=et[v[4]]*tk0[v[4]]*ro[v[4]]*(eps[11]*GXKOEF+eps[12]*GYKOEF);
			ival2=et[v[4]]*tk0[v[4]]*(eps[11]*dpdx+eps[12]*dpdy);
			fprintf(fp_log,"%ld %ld   %e %e %e %e   %e %e %e\n",m1,m2,e,n,dpdx,dpdy,et[v[4]],ival1,ival2); getchar();
			*/
		}
		/* Simplified adiabatic term alphV*T/Cp*(gY*vY+gX*vX) */
		if(adiabyn==1)
		{
			ival+=et[v[4]]*tk0[v[4]]/cp[v[4]]*(eps[11]*GXKOEF+eps[12]*GYKOEF);
		}

		/* Viscouse friction dH/dt=NU*EPSik*dVi/dXk */
		/* Viscouse friction dH/dt=2NU*EPSik*EPSik=SIGik*EPSik, where SIGik=2NU*EPSik, EPSik=1/2(dVxi/dXk+dVxk/dXi) */
		/* Nu values Use */
		/* Check Vx, Vy Boundary conditions  around T Node */
		// Set to 255 in mode.t3c to avoid shear heating near the boundary conditions applied before node 255
		if(frictyn && (m1-1)>=frictyn)
			if (!bondm[((m1-frictyn)*ynumy+m2)*3+1] && !bondm[(m1*ynumy+m2+frictyn)*3+1])
				if (tk0[v[4]]>0.0)
		{
			/* EPSxx*SIGxx  EPSyy*SIGyy  EPSxy*SIGxy interpolated values  */
			epsval=2.0*eps[13]+2.0*eps[14];
			/**/
			ival+=epsval/cp[v[4]]/ro[v[4]];
			/**/
			/*
			fprintf(fp_log,"A %ld %ld   %e\n",m1,m2,epsval); getchar();
			fprintf(fp_log,"%ld %ld %ld %e %e %e %e",m1,m2,v[4],tk0[v[4]],kt[v[4]],cp[v[4]],ro[v[4]]); getchar();
			*/
		}

		/* Save Right Part for next calculations */
		sol1[v[4]]=ival;

		/* 0       3       6 */
		/*                   */
		/*        Qy3        */
		/*                   */
		/* 1  Qx1  4  Qx4  7 */
		/*       <T4>        */
		/*                   */
		/*        Qy4        */
		/*                   */
		/* 2       5       8 */
		/*  1/ro/Cp*[dQx/dX + dQy/dY] */
		ival+=((qxcalc(m1,m2,0)-qxcalc(m1-1,m2,0))/xkf+(qycalc(m1,m2,0)-qycalc(m1,m2-1,0))/ykf)/cp[v[4]]/ro[v[4]];
		/**/
		/* Rate of heating return */
		return ival;

		/* Left Side Calculation ---------------------- */
		case 1:
		/* Implicit form */
		/* dT/dt - 1/Cp/ro*[dQx/dX + dQy/dY] =  1/Cp/ro*(dH/dt] */
		/**/
		/* Initial position Number */
		wn[0]=1;
		/**/
		/* Right part  1/Cp/ro*(dH/dt) */
		wi[0]=sol1[v[4]]+tk0[v[4]]/timestep;
		/**/
		/* Add Left part */
		/* dT/dt */
		wn[1]=v[4];
		wi[1]=1.0/timestep;
		/* - 1/Cp/ro*[dQx/dX] */
		qxcalc(m1-1,m2,1.0/xkf/cp[v[4]]/ro[v[4]]);
		qxcalc(m1,m2,-1.0/xkf/cp[v[4]]/ro[v[4]]);
		/* - 1/Cp/ro*[dQy/dY] */
		qycalc(m1,m2-1,1.0/ykf/cp[v[4]]/ro[v[4]]);
		qycalc(m1,m2,-1.0/ykf/cp[v[4]]/ro[v[4]]);
		/**/
		return 0;

		/* Error Calculation ---------------------- */
		case 2:
		/**/
		/* Left part Add */
		/* dT/dt */
		ival=(tk[v[4]]-tk0[v[4]])/timestep;
		/**/
		/*  -1/ro/Cp*[dQx/dX + dQy/dY] */
		ival-=((qxcalc(m1,m2,0)-qxcalc(m1-1,m2,0))/xkf+(qycalc(m1,m2,0)-qycalc(m1,m2-1,0))/ykf)/cp[v[4]]/ro[v[4]];
		/**/
		/* Right part Add */
		ival-=sol1[v[4]];
		return ival;
	}
	return 0;
}
/* Left side, Right side or Err of Thermal Conductivity equation calc */



/* Coefficients or value for Qx  Equation */ 
/* Qx=k(dT/dX) */
double qxcalc(long int m1, long int m2, double ynval)
	/* m1,m2 - node X,Y number */
	/* ynval - Val Qx Calc Y(0)/N (koefficient) */
{
	/* Qx horizontal position */
	double xi=(gx[m1]+gx[m1+1])/2.0,leftqx=0,kteff;
	long int m1min,m1max,m3;
	int n1,n;

	/* Effective heat conductivity calc */
	kteff=(kt[m1*ynumy+m2]+kt[(m1+1)*ynumy+m2])/2.0;

	/* dT/dX */
	/* Calc, Check Fd limits */
	m1min=m1-heatfd;
	if(m1min<0) m1min=0;
	m1max=m1+1+heatfd;
	if(m1max>xnumx-1) m1max=xnumx-1;

	/* Load distances to xn[] */
	for (m3=m1min;m3<=m1max;m3++)
	{
		xn[m3-m1min]=gx[m3];
	}

	/* Calc maximal position in xn[] */
	n=(int)(m1max-m1min);

	/* Calc T coefficients for Qx */
	fdweight(n,1,xi);

	/* Return val for LSQ err ----------------------------*/
	if(ynval==0)
	{
		/* RIGHT part of Qx */
		/* Qx=k(dT/dX) */
		/* k(dT/dX)=0 */
		/* Add T with koefficients */
		for (m3=m1min;m3<=m1max;m3++)
		{
			leftqx+=tk[m3*ynumy+m2]*cn[m3-m1min][1]*kteff;
		}

		return leftqx;
	}

	/* Add Coefficients for left parts of Qx ----------------*/
	/* k(dT/dX)=0 */
	/* Add T with koefficients */
	for (m3=m1min;m3<=m1max;m3++)

	{
		wn[wn[0]+1+m3-m1min]=m3*ynumy+m2;
		wi[wn[0]+1+m3-m1min]=ynval*cn[m3-m1min][1]*kteff;
	}

	/* Add total Num of lines -------------------------------------- */
	wn[0]+=n+1;

	return 0;
}
/* Coefficients or value for Qx  Equation */ 





/* Coefficients or value for Qy  Equation */ 
/* Qy=k(dT/dY) */
double qycalc(long int m1, long int m2, double ynval)
	/* m1,m2 - node X,Y number */
	/* ynval - Val Qy Calc Y(0)/N(koefficient) */
{
	/* Qx horizontal position */
	double xi=(gy[m2]+gy[m2+1])/2.0,leftqy=0,kteff;
	long int m2min,m2max,m3;
	int n1,n;

	/* Effective heat conductivity calc */
	kteff=(kt[m1*ynumy+m2]+kt[m1*ynumy+m2+1])/2.0;

	/* dT/dY */
	/* Calc, Check Fd limits */
	m2min=m2-heatfd;
	if(m2min<0) m2min=0;
	m2max=m2+1+heatfd;
	if(m2max>ynumy-1) m2max=ynumy-1;

	/* Load distances to xn[] */
	for (m3=m2min;m3<=m2max;m3++)
	{
		xn[m3-m2min]=gy[m3];
	}

	/* Calc maximal position in xn[] */
	n=(int)(m2max-m2min);

	/* Calc T coefficients for Qy */
	fdweight(n,1,xi);

	/* Return val for LSQ err ----------------------------*/
	if(ynval==0)
	{
		/* RIGHT part of Qy */
		/* Qy=k(dT/dY) */
		/*  k(dT/dY)=0 */
		/* Add T with koefficients */
		for (m3=m2min;m3<=m2max;m3++)
		{
			leftqy+=tk[m1*ynumy+m3]*cn[m3-m2min][1]*kteff;
		}

		return leftqy;
	}

	/* Add Coefficients for left parts of Qx ----------------*/
	/* -k(dT/dY)=0 */
	/* Add T with koefficients */
	for (m3=m2min;m3<=m2max;m3++)

	{
		wn[wn[0]+1+m3-m2min]=m1*ynumy+m3;
		wi[wn[0]+1+m3-m2min]=ynval*cn[m3-m2min][1]*kteff;
	}

	/* Add total Num of lines -------------------------------------- */
	wn[0]+=n+1;

	return 0;
}
/* Coefficients or value for Qy  Equation */ 




/* Left side or Err for T Boundary Condition Equation */ 
/* T=CONST+KOEF*Tn */
double tbonderr(long int mcmax, int ynerr)
	/* mcmax - numer of cur T in sol[] */
	/* ynerr - Err Calc Y(1)/N(0) */
{
	/* Val Buffer */
	double leftt=0;
	/* Boundary condition number for T */
	long int mcmax1=mcmax+nodenum3;
	int n1;

	/* Error Calc */
	if (ynerr)
	{
		/* T Simple boundary conditions */
		/* Add Const */
		leftt=x[mcmax]-bondv[bondm[mcmax1]][0];
		/* Add Koef */
		if(bondn[bondm[mcmax1]][0]) leftt-=bondv[bondm[mcmax1]][1]*x[bondn[bondm[mcmax1]][0]-1-nodenum3];
		/* Add Koef1 */
		if(bondn[bondm[mcmax1]][1]) leftt-=bondv[bondm[mcmax1]][2]*x[bondn[bondm[mcmax1]][1]-1-nodenum3];
		/* Add Koef2 */
		if(bondn[bondm[mcmax1]][2]) leftt-=bondv[bondm[mcmax1]][3]*x[bondn[bondm[mcmax1]][2]-1-nodenum3];
		return leftt;
	}

	/* Add T CONST */
	wn[0]=1;
	wi[0]=bondv[bondm[mcmax1]][0];
	wn[1]=mcmax;
	wi[1]=1.0;
	/* Add T PAR1,PAR2,PAR3 */
	for (n1=0;n1<3;n1++)
	{
		if(bondn[bondm[mcmax1]][n1]) 
		{
			wn[0]+=1;
			wn[wn[0]]=bondn[bondm[mcmax1]][n1]-1-nodenum3;
			wi[wn[0]]=-bondv[bondm[mcmax1]][n1+1];
		}
	}

	return 0;
}
/* Left side or Err for T Boundary Condition Equation */ 



/* tk[] Correct for boundary conditions */
void tkrecalc()
{
	/* Counters */
	long int m1,m2,mcmax,mcmax1;
	int n1;

	/* Save Old temperatures  for Bondary solitions */
	pos0cur=0;
	/* Save TK */
	for (m1=0;m1<nodenum;m1++)
	{
		tk1[m1]=tk[m1];
		sol0[m1]=0;
		num0[m1]=0;
	}

	/* Add Matrix by Boundary Equations */
	if (printmod) fprintf(fp_log,"Adding Matrix By TK Boundaries ...");
	fflush(fp_log);
	/* TK - Node Cycle */
	for (m1=0;m1<xnumx;m1++)
		for (m2=0;m2<ynumy;m2++)
	{
		/* Cur Line Num in sol0[] */
		mcmax=m1*ynumy+m2;
		mcmax1=mcmax+nodenum3;
		if(!bondm[mcmax1])
		{
			/* Thermal conductivity Equation koef calc add */
			wn[0]=1;
			wi[0]=tk1[mcmax];
			wn[1]=mcmax;
			wi[1]=1.0;

			/* Add matrix */
			// ATAT TARAS: do we need this gausmat3 (below 3x, are only ones)? could it be done with gausmat4 instead?
			gausmat3(3,mcmax,0);
		}
		else	
		{
			/* Add vX Simple Boundary */
			tbonderr(mcmax,0);
			/**/
			/* Add matrix */
			gausmat3(3,mcmax,0);
		}
	}

	/* Solve Matrix */
	mcmax=nodenum-1;
	gausmat3(0,mcmax,0);
	/* End Solve Matrix */

	/* TK reload */
	for (m1=0;m1<nodenum;m1++)
	{
		tk[m1]=sol0[m1];
	}
	if (printmod) fprintf(fp_log,"OK!\n");
	fflush(fp_log);
	/* End Boundary Solution Cycle ------------------------------------------ */
}
/* End tk[] Correct for boundary conditions */


