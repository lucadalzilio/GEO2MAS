//** 2 functions called from i2.c **//

/* Move markers by using simple Runge-Kutta method */
void movemarkomp()
{
	/* Vx, Vy buffer */
	double dvxdx,dvxdy,dvydx,dvydy,celdx,celdy,vx0,vx1,vx2,vx3,vx4,vy0,vy1,vy2,vy3,vy4,ee0,ee1,ee2,ee3,ee4,sp0,sp1,sp2,sp3,sp4,pr0,pr1,pr2,pr3,pr4;
	/* Water */
	double vxwater,vywater;
	long int mm1,marknum1,m10,m20,m30,m1,m2,m3;
	/* Erosion-Sedimentation Y/N */
	int n1;
	int mm2;
	/* Nonstabilyty for immobile markers */
	double xnonstab=0.50,ynonstab=0.60;
	double dpdx,dpdy,e,n,vxkoef,vykoef,dx,dy;
	double start;

	/* Hydration front progress */
#if setup>9
	start=omp_get_wtime();	
	/* dehydration */ 
	if(vyfluid!=0 && timesum>1e+11) hydration2omp();
	fprintf(fp_log,"\n  Time taken for hydration = %e s \n",omp_get_wtime()-start);
	
	
#endif

	/* Save number of markers */
	marknum1=marknum;

	/* Surface changes */
#if setup>9
	start=omp_get_wtime();	
	if (timestep && erosmod) erosion();
	fprintf(fp_log,"\n  Time taken for erosion = %e s \n",omp_get_wtime()-start);
#endif
	
	/* Move markers */	
#pragma omp parallel for shared(markx,marky,markt,markim,markk,markxx,markxy,markd,markv,markp,markexx,markexy,markw,marke,marknum1,follow,nm,ystpy,m10_hr,m11_hr,marknum,outgrid,eroslev,vyfluid,vymelt,GXKOEF,GYKOEF,markmod,timestep,zdeep,tdeep,markht,xsize,ysize,xnumx,ynumy,gx,gy,pr,esp,exx,exy,vx,vy,markcp,markkt,markkf,markkp,markro,markbb,markaa,marknu,markn0,markn1,markll,marka0,marka1,markb0,markb1,markdh,markdv,markss,markmm,marks1,start_cond,stoksmod) \
	private(mm1,mm2,m10,m20,m1,m2,m3,n1,vxwater,vywater,e,n,dpdx,dpdy,a,b,vxkoef,vykoef,vx0,vy0,sp0,ee0,pr0,vx1,vy1,sp1,ee1,pr1,vx2,vy2,sp2,ee2,pr2,vx3,vy3,sp3,ee3,pr3,vx4,vy4,sp4,ee4,pr4,dx,dy) \
	schedule(runtime)
	for (mm1=0;mm1<marknum;mm1++)
	{
		/* Marker type */
		mm2=(int)markt[mm1]; if (mm2>=100) mm2-=100;
		if( ((markx[mm1]>=0 && marky[mm1]>=0 && (markx[mm1])<=xsize && (marky[mm1])<=ysize) || outgrid!=1) && !markim[mm2] )
		{
			// Search marker location within nodal grid
			m10=m1serch(markx[mm1]);
			m20=m2serch(marky[mm1]);
			/**/
			/* Erosion-Sedimentation */
			if((marky[mm1])<=eroslev) n1=1; else n1=0;

			/* Water marker move */
			vxwater=vywater=0;
			if(markt[mm1]>=50 && markt[mm1]<100) 
			{
				/* Water velocity */
				vywater=vyfluid; if(markd[mm1]>1100.0) vywater=vymelt;
				/* Fluid in rock */
				if(vyfluid>0 && (markk[mm1]==0 || markk[mm1]>298.0)) 
				{
					/* Horizontal,Vertical P-cell index */
					m1=m10; if(markx[mm1]>(gx[m1]+gx[m1+1])/2.0) m1+=1;
					if(m1<1) m1=1; if(m1>xnumx-2) m1=xnumx-2;
					m2=m20; if(marky[mm1]>(gy[m2]+gy[m2+1])/2.0) m2+=1;
					if(m2<1) m2=1; if(m2>ynumy-2) m2=ynumy-2;
					/* Pressure gradients */
					e=(markx[mm1]-(gx[m1-1]+gx[m1])/2.0)/((gx[m1+1]-gx[m1-1])/2.0);
					n=(marky[mm1]-(gy[m2-1]+gy[m2])/2.0)/((gy[m2+1]-gy[m2-1])/2.0);
					m3=m1*ynumy+m2;
					dpdx=2.0*((1.0-n)*(pr[m3+ynumy]-pr[m3])+n*(pr[m3+ynumy+1]-pr[m3+1]))/(gx[m1+1]-gx[m1-1]);
					dpdy=2.0*((1.0-e)*(pr[m3+1]-pr[m3])+e*(pr[m3+ynumy+1]-pr[m3+ynumy]))/(gy[m2+1]-gy[m2-1]);
					/* Recalc velocity koefficients */
					vxkoef=(1000.0*GXKOEF-dpdx)/(2300.0*9.81);
					vykoef=(1000.0*GYKOEF-dpdy)/(2300.0*9.81);

					if(vxkoef>2.0) vxkoef=2.0; if(vxkoef<-2.0) vxkoef=-2.0;
					if(vykoef>2.0) vykoef=2.0; if(vykoef<-2.0) vykoef=-2.0;
					/* Recalc velocity */
					vxwater=vywater*vxkoef;
					vywater*=vykoef; 
				}
				else
					/* Fluid in water */
				{
					vxwater=0;
					vywater=-ABSV(vywater);
				}
			}
	
			/* Motion Calc ///////////////////////////////// */
	
			/* Vx, Vy, EpsII Simple calc */
			if(markmod==1)
			{
				/* Interpolate velocity, pressure?, EE(eii), and ESP (spin) */
				// These marker values are not stored, they are used in this routine to move each marker in private, thats it.	
				allinteriomp(markx[mm1],marky[mm1],m10,m20,&vx0,&vy0,&pr0,&sp0,&ee0); 
				vx0+=vxwater; vy0+=vywater;
				/**/
				/* fprintf(fp_log,"SIMPLE %ld %d %e %e   %e %e %e",mm1,markt[mm1],markx[mm1],marky[mm1],vx0,vy0,sp0); getchar(); */
			}
			/* Vx, Vy, EpsII 4 Runge-Kutta koef calc */
			else
			{
				allinteriomp(markx[mm1],marky[mm1],m10,m20,&vx1,&vy1,&pr1,&sp1,&ee1);
				vx1+=vxwater; vy1+=vywater;
				/**/
				//fprintf(fltest,"RK4   %ld %d %e %e   %e %e %e %e \n",mm1,markt[mm1],markx[mm1],marky[mm1],vx1,vy1,sp1,ee1);
				/**/
				allinteriomp(markx[mm1]+vx1*timestep/2.0,marky[mm1]+vy1*timestep/2.0,m10,m20,&vx2,&vy2,&pr2,&sp2,&ee2);
				vx2+=vxwater; vy2+=vywater;
				/**/
				allinteriomp(markx[mm1]+vx2*timestep/2.0,marky[mm1]+vy2*timestep/2.0,m10,m20,&vx3,&vy3,&pr3,&sp3,&ee3);
				vx3+=vxwater; vy3+=vywater;
				/**/
				allinteriomp(markx[mm1]+vx3*timestep,marky[mm1]+vy3*timestep,m10,m20,&vx4,&vy4,&pr4,&sp4,&ee4);
				vx4+=vxwater; vy4+=vywater;
				/**/
				/* Vx,Vy, EpsXX, EpsYY, EpsXY calc after Runge-Kutta */
				vx0=(vx1+2.0*vx2+2.0*vx3+vx4)/6.0;
				vy0=(vy1+2.0*vy2+2.0*vy3+vy4)/6.0;
				if(markmod==2)
				{
					sp0=(sp1+2.0*sp2+2.0*sp3+sp4)/6.0;
					ee0=(ee1+2.0*ee2+2.0*ee3+ee4)/6.0;
				}
				else
				{
					sp0=sp1;
					ee0=ee1;
				}
			}
	
			/* Orthogonal motion only */
			if (outgrid==2)
			{
				if(markx[mm1]<0 || (markx[mm1])>xsize) vy0=0;		
				if(marky[mm1]<0 || (marky[mm1])>ysize) vx0=0;		
			}

			/* Normal markers */
			if(markt[mm1]<100)
			{
				/* Markers coming from below the model */
				if(marky[mm1]>zdeep && markk[mm1]<tdeep) markk[mm1]=tdeep;
				// If you do not want to apply a large temperature lower boundary condition use:
				//if(marky[mm1]>zdeep && vy0<0 && markk[mm1]<tdeep) markk[mm1]=tdeep;
				
				/* Normal markers */
				/* X,Y calc after Runge-Kutta */
				markx[mm1]+=(timestep*vx0);
				marky[mm1]+=(timestep*vy0);
				if(marke[mm1]>0)
				{
					marke[mm1]+=(timestep*ee0);
				}
				sp0*=timestep;
				
				/* Turcotte & Schubert, 1995 rotation formula */
				if(stoksmod==1)
				{
					sp1=markxx[mm1]*cos(sp0)*cos(sp0)-markxx[mm1]*sin(sp0)*sin(sp0)+markxy[mm1]*sin(2.0*sp0);
					sp3=0.5*(-markxx[mm1]-markxx[mm1])*sin(2.0*sp0)+markxy[mm1]*cos(2.0*sp0);
					markxx[mm1]=sp1;
					markxy[mm1]=sp3;
				}
				
				/* Jaumann corrotation formula */
				if(stoksmod==2)
				{
					sp1=markxx[mm1]+markxy[mm1]*2.0*sp0;
					sp3=markxy[mm1]+0.5*(-markxx[mm1]-markxx[mm1])*2.0*sp0;
					markxx[mm1]=sp1;
					markxy[mm1]=sp3;
				}

				/* Out of grid marker reset */
				if(markx[mm1]<0 || marky[mm1]<0 || (markx[mm1])>xsize || (marky[mm1])>ysize) 
				{
					markk[mm1]=0;
					markd[mm1]=-1.0;
					markw[mm1]=-1.0;
					marke[mm1]=0;
				}
			}
			/* Immobile markers */
			else
			{
				/* X,Y calc after Runge-Kutta */
				// Which velocity is used here, if were before located outside of grid ...
				markx[mm1]+=(timestep*vx0);
				marky[mm1]+=(timestep*vy0);
				
				/* Check new position, add marker at end (marknum1) */
				// Immobile markers that now enter grid
				if(markx[mm1]>=0 && marky[mm1]>=0 && markx[mm1]<=xsize && marky[mm1]<=ysize)
				{
#pragma omp critical(newmark)
					{
#pragma omp flush(marknum1)
						/* Type save */
						markt[marknum1]=markt[mm1]-100;
						/* X,Y calc after Runge-Kutta */
						// Give marker new location (within grid)
						markx[marknum1]=markx[mm1];
						marky[marknum1]=marky[mm1];
						/* Temperature Reset */
						markk[marknum1]=0;
						markd[marknum1]=-1.0;
						markv[marknum1]=0;
						/* Strain Reset */
						marke[marknum1]=0;
						/* Stress Reset */
						markxx[marknum1]=0;
						markxy[marknum1]=0;
						/* Pressure Reset */
						markp[marknum1]=0;
						/* Strain rate Reset */
						markexx[marknum1]=0;
						markexy[marknum1]=0;
						/* Add aditional markers counter */
						marknum1++;

						/* X,Y reset for immobile marker */
						markx[mm1]=markk[mm1];
						marky[mm1]=markv[mm1];

						// If new marker is interesting for picking algorithm, flag to follow
						// Note is hard-coded in i2.c as well. Only here excluded fluid markers, since immobile can not become fluid
#if setup>9
						if (start_cond==1 && marky[marknum1]<85e3 && markx[marknum1]>gx[m10_hr] && markx[marknum1]<gx[m11_hr] && markt[marknum1]>1 && markt[marknum1]<50)
						{ 
							follow[marknum1]=1;
							// #pragma omp flush(nm)
							nm++;
						}
#endif
					}
				}
				/* Check,Reset old position */
				// Use markk and v as dummy from above, so dx and/or dy are 0 if marker is newly added
				dx=markx[mm1]-markk[mm1];
				dy=marky[mm1]-markv[mm1];
				dy=pow(dx*dx+dy*dy,0.5);
				/*
				if(dy>ystpy || (marky[mm1]<0 && vy0<0) || (marky[mm1]>ysize && vy0>0) || (markx[mm1]<0 && vx0<0) || (markx[mm1]>xsize && vx0>0))
				*/
				// If moved by more than one cell, reset to old position ?
				if(dy>ystpy)
				{
					/* X,Y reset for immobile marker */
					markx[mm1]=markk[mm1];
					marky[mm1]=markv[mm1];
				}
			}
			/* End Motion Calc ///////////////////////////////// */
		}
	}
	// End omp-section move markers

	/* Mark num */
	if(marknum1>MAXMRK) {fprintf(fp_log,"Space out in markx[]"); fflush(fp_log); exit(0);}

	/* Reset aditional markers */
	mm1=0;
	while(marknum1>marknum && mm1<marknum)
	{
		/* Reload marker */
		if((markx[mm1]<0 || marky[mm1]<0 || (markx[mm1])>xsize || (marky[mm1])>ysize) && markt[mm1]<100) 
		{
			/* Decrease aditional markers counter */
			marknum1--;
			/* Type save */
			markt[mm1]=markt[marknum1];
			/* Temperature Reset */
			markk[mm1]=0;
			markd[mm1]=-1.0;
			/* Strain Reset */
			marke[mm1]=0;
			/* Stress Reset */
			markxx[mm1]=0;
			markxy[mm1]=0;
			/* Pressure Reset */
			markp[mm1]=0;
			/* Strain rate Reset */
			markexx[mm1]=0;
			markexy[mm1]=0;
			/* X,Y reload  */
			markx[mm1]=markx[marknum1];
			marky[mm1]=marky[marknum1];
		}
		/* Increase markers counter */
		mm1++;
	}
	fprintf(fp_log,"\n Number of markers: OLD = %ld NEW = %ld \n",marknum,marknum1);fflush(fp_log);
	
	/* Set new marker number */
	marknum=marknum1;

	/* Incr cycle of sedimentation */
	sedimnum++;
}
/* End OMP move markers by using Simple/Runge-Kutta method */



/* ro[],nu[] recalc after marker positions */
void ronurecalcomp()
{
	/* Counters */
	long int m1,m2,m3,m10,m20; 
	int mm2,yn,mm3,n1,n2,ncount=0,nt,tid;
	long int mm1;
	double dx,dy,swt,swt1,celdx,celdy;
	double wro,mnu,mgg,maa,mdro,msxxe,msxye,mexxe,mexye,mro,mcp,mkt,mht,mbb,mdi0,mdi1,mwa,dmwa,mxmelt,mhlatent;
	double Mgg,Mro,Mwa,Mcp,Mbb,Maa,Mdhh,Mkt;
	// Here in this loop epsin is 2nd invariant of visco-plastic strainrate
	double sigin,epsin;
	/* TD Database variables,  dTK,dPB - TK, PB step for tabulation in TD database */
	double H0,H1,H2,H3,R0,R1,R2,R3,G0,G1,G2,G3,W0,W1,W2,W3,dTK=20.0,dPB=1000.0,n,e;
	/* Phase transition variables */
	double p_pl_out,p_ga_in,rokf,p_sp_in,p_ol_out,p_pv_in,p_sp_out,p_st_in;
	/* RO, NU equations var */
	double mpb=1.0,mtk=300.0,numax=0,numin=0;
	double start,xwall,b1,b2,slope_wall,gelbeg;
	
	start=omp_get_wtime();
	
	Mgg=Mro=Mwa=Mcp=Mbb=Maa=Mdhh=Mkt=0;
	
	if (printmod) fprintf(fp_log,"\n Number of nodes = %ld  Number of markers = %ld \n",nodenum,marknum);
	fflush(fp_log);	
	
#pragma omp parallel
	{nt=omp_get_num_threads();}
	
	/* Layering on sediments */
	m1=(long int)(sedimnum/sedimcyc);
	m2=((long int)(m1/2))*2;
	if(m2==m1) yn=3; else yn=4;
	
	/* ADD MARKERS TO THE v-CELLS ========================== */
	/* Clear ro[],nu[] wt */
	for (m1=0;m1<nodenum;m1++)
	{
		ro0[m1]=0;
		et0[m1]=0;
		nu0[m1]=0;
		nd0[m1]=0;
		gg0[m1]=0;
		gd0[m1]=0;
		sxxe0[m1]=0;
		sppe0[m1]=0;
		sbritn0[m1]=0;  	// yield stress
		sxye0[m1]=0;
		exxe0[m1]=0;
		exye0[m1]=0;
		dro0[m1]=0;
		drp0[m1]=0;
		cp0[m1]=0;
		kt0[m1]=0;
		ht0[m1]=0;
		tk0[m1]=0;
		mrx0[m1]=0;
		mry0[m1]=0;
		mvx0[m1]=0;
		mvy0[m1]=0;
		sol0[m1]=0;
		sol0[nodenum+m1]=0;
		sol0[nodenum2+m1]=0;
		sol1[m1]=0;
		sol1[nodenum+m1]=0;
		sol1[nodenum2+m1]=0;
	}
		
#if setup>9
	/* (1) Erosion-sedimentation and melting account for all markers */
#pragma omp parallel for shared(markx,marky,markk,markt,marke,markd,marknum,waterlev,erosmod,deserp,dyserp,xsize,ysize,gx,gy,xnumx,ynumy,ep,pr,timesum,res_high) \
	private(mm1,mm2,m3,mtk,mpb,m10,m20,mxmelt,mhlatent) \
	firstprivate(yn) \
	schedule(runtime)
	for (mm1=0;mm1<marknum;mm1++)
	{
		/* Check markers out of grid */
		if(markx[mm1]>0 && marky[mm1]>0 && (markx[mm1])<xsize && (marky[mm1])<ysize && markk[mm1]>0 && markt[mm1]<50)
		{

			/* Up Left Node X,Y Num */
			m10=m1serch(markx[mm1]);
			m20=m2serch(marky[mm1]);
			m3=m10*ynumy+m20;	
			
			mm2=(int)markt[mm1];
			
			/* Erosion/sedimentation account */
			if(erosmod) erosmarkomp(mm1,yn,m10,markx[mm1],marky[mm1],markt,marke,markd);
			
			/* Water/Air account */
			if(markt[mm1]<2)
			{
				/* Change marker type */
				if((marky[mm1])>waterlev) markt[mm1]=1; else markt[mm1]=0;
			}
			
			/* P, T parameters calc */
			// 1e-5 since convert to bars for melting look-up?
			mpb=1e-5*allinterpomp(markx[mm1],marky[mm1],m10,m20);
			mtk=markk[mm1];
		
			// Remove initial weak zone for subduction initiation after X My
			if(timesum>(20.0e6*3.15576e+7) && markt[mm1]==12)
				{markt[mm1]=9;}

			/* Serpentinization of brittle mantle faults at sub-surface */
			if((markt[mm1]==9 || markt[mm1]==9 || markt[mm1]==9) && marke[mm1]>deserp && marky[mm1]<dyserp) 
			{
				/* Mantle to Antigorite transformation */
				markt[mm1]=13; 
				markd[mm1]=-1.0; 
			}
			
			/* Mantle to Antigorite transformation */
			antigoromp(mtk,mpb,markx[mm1],marky[mm1],mm1,m10,markt);
			
			/* Rocks to rock+melt transformation */
			// Note markt passes in address of first element of array to function and allows for modification there
			if (meltmod) meltingomp(mtk,mpb,mm1,mm2,markt,marke,&mxmelt,&mhlatent);
		}
	}
	// End OMP section erosion-sedimentation
	
	// Open file for storing marker interface data
	if (n0==1)
	{
		flfric = fopen(fileTxtOutputFric,"w");
		fprintf(flfric," rocktype markx marky markpressure sbrit markxx markxy \n");
	}
#endif

	if (printmod==10000) fprintf(fp_log,"\n Time taken for erosmark/antigor/melting in ronurecalc = %e s \n",omp_get_wtime()-start);
	start=omp_get_wtime();
		
	/* (2) Add ro[] nu[] etc. using selected markers */
#pragma omp parallel shared(marknum,gridmod,markx,marky,markk,markt,erosmod,sedilev,markw,markd,eroslev, \
	waterlev,zdeep,densimod,markht,marke,markexx,markexy,markxx,markxy,markp, \
	nodenum,nodenum2,markvx,markvy,markv,markwa,markim,xsize,ysize,gx,gy,xnumx,ynumy,ep,pr,vx,vy,markf0,markf1,markbb,markaa,markro,markgg, \
	markn0,markn1,marks0,marks1,marknu,markcp,markkt,markkf,markkp,nubeg,nuend,strmin,strmax,\
	exy,esp,exx,pbmin,pbstp,pbnum,tkmin,tkstp,tknum,pbmin1,pbstp1,pbnum1,tkmin1,tkstp1,tknum1,timesum,td, \
	zmpor,tkpor,markll,hidrl,hidry,lambfld,marka0,markb0,marke1,marka1,markb1,marke0,msbrit,msii_old, \
	tk_updipsez0,tk_updipsez1,mgamma_vw,mgamma_vs,mvc_vs,mvc_vw,mus_vs,markdh,markdv,markss,markmm,timestepe,cyc0max,start_cond,veldepfric,stoksmod,res_high) \
	private(mm1,mm2,tid,m10,m20,m1,m3,mpb,mtk,mro,mbb,maa,mcp,mkt,mht,mnu,mgg,mxmelt,mhlatent,mwa,wro,dmwa,mdi0,mdi1,mdro, \
	celdx,celdy,swt,swt1,dx,dy,msxxe,msxye,mexxe,mexye,sigin,epsin,Mgg,Mro,Mwa,Mcp,Mbb,Maa,Mdhh,Mkt,p_pl_out,p_ga_in,rokf,p_sp_in,p_ol_out,p_pv_in,p_sp_out,p_st_in) \
	firstprivate(yn)
	{				
		// Initialize temporarily interpolation arrays (capitilized) to zero (inside pragma, so is private)
		double* Nu0 	= (double*) calloc(nodenum,sizeof(double));
		double* Nd0 	= (double*) calloc(nodenum,sizeof(double));
		double* Gg0 	= (double*) calloc(nodenum,sizeof(double));	
		double* Gd0 	= (double*) calloc(nodenum,sizeof(double));
		double* Ro0 	= (double*) calloc(nodenum,sizeof(double));	
		double* Sxxe0 	= (double*) calloc(nodenum,sizeof(double));
		double* Sppe0 	= (double*) calloc(nodenum,sizeof(double));
		double* Sbritn0 = (double*) calloc(nodenum,sizeof(double));	
		double* Sxye0 	= (double*) calloc(nodenum,sizeof(double));
		double* Exxe0 	= (double*) calloc(nodenum,sizeof(double));
		double* Exye0 	= (double*) calloc(nodenum,sizeof(double));
		double* Et0 	= (double*) calloc(nodenum,sizeof(double));
		double* Dro0 	= (double*) calloc(nodenum,sizeof(double));	
		double* Drp0 	= (double*) calloc(nodenum,sizeof(double));
		double* Cp0 	= (double*) calloc(nodenum,sizeof(double));
		double* Kt0 	= (double*) calloc(nodenum,sizeof(double));
		double* Ht0 	= (double*) calloc(nodenum,sizeof(double));	
		double* Tk 	= (double*) calloc(nodenum,sizeof(double));
		double* Mrx0 	= (double*) calloc(nodenum,sizeof(double));
		double* Mry0 	= (double*) calloc(nodenum,sizeof(double));
		double* Mvx0 	= (double*) calloc(nodenum,sizeof(double));	
		double* Mvy0 	= (double*) calloc(nodenum,sizeof(double));
		double* Sol0 	= (double*) calloc(nodenum*3,sizeof(double));
		double* Sol1 	= (double*) calloc(nodenum*3,sizeof(double));
	
#pragma omp for \
		schedule(runtime)
		for (mm1=0;mm1<marknum;mm1+=gridmod)
		{	
			/* Check markers out of grid */
			if(markx[mm1]>0 && marky[mm1]>0 && (markx[mm1])<xsize && (marky[mm1])<ysize && markk[mm1]>0 && markt[mm1]<50)
			{
				tid=omp_get_thread_num();
			
				m10=m1serch(markx[mm1]);
				m20=m2serch(marky[mm1]);
			
				/* Marker type */
				mm2=(int)markt[mm1]; if (mm2>=100) mm2-=100;
			
#if setup>9
				/* 1a. --- Remove water, rocks --- */
				if(erosmod==0)
				{
					if(marky[mm1]>sedilev && mm2<2) 
					{
						mm2=yn; markt[mm1]=yn;
						markw[mm1]=0;
						markd[mm1]=-1.0;
					}
					if(marky[mm1]<eroslev && mm2>1) 
					{
						if((marky[mm1])>waterlev) markt[mm1]=1; else markt[mm1]=0;
						mm2=markt[mm1];
						markw[mm1]=0;
						markd[mm1]=-1.0;
					}
				}

				/* 1b. Remove Plumes */
				if(marky[mm1]>zdeep && mm2!=10) 
				{
					mm2=10; markt[mm1]=10;
					markw[mm1]=0;
					markd[mm1]=-1.0;
				}
#endif

				/* P, T parameters calc */
				mpb=1e-5*allinterpomp(markx[mm1],marky[mm1],m10,m20);
				mtk=(markk[mm1]);
			
				/* Reset water/air temperature */
				// if (mm2<2) mtk=markk[mm1]=273.0;

				/* 2.-3. --- Calculate density --- */
				// & just points to address of normally defined variables at start of this routine; more efficient for passing! this way variable here can be changed inside subroutine
#if setup>9
				dencalcomp(mtk,mpb,markx[mm1],marky[mm1],mm2,&mro,&mbb,&maa);
				mcp=markcp[mm2];
				mkt=(markkt[mm2]+markkf[mm2]/(mtk+77.0))*exp(markkp[mm2]*mpb);

				/* ========================= */
				/* Mantle phase transitions  */
				/* ========================= */
				if (densimod==1)
				{
			        /*
				if(mm2>=9 && mm2<=14 && markex[mm1]>0) mro*=1.0-0.04*markex[mm1];
		    	        */
			  	    /* Eclogitization, St, Pv transitions in oceanic crust */
			        if(mm2==7 || mm2==8 )
			        {
				    /* Eclogitization Ito and Kennedy, 1971 */
			        /*basalt=>garnet granulite (Ga-In) transition*/
			        p_ga_in=-9222.0+mtk*14.0;
			        /*Not to have granulites at pressure lower than 2 kbar*/
			        if(p_ga_in<2000.0) p_ga_in=2000.0;
			        /*garnet granulite=>eclogite (Pl-Out) transition*/
			        p_pl_out=-1460.0+mtk*20.0;
			        /*Not to have eclogites at pressure lower than 12 kbar*/
			        if(p_pl_out<12000.0) p_pl_out=12000.0;
			        if(mpb>p_ga_in)
			                {
					rokf=0;
					if(mtk>teclmin)
						{
						if(mtk>teclmax)
			       				{
			       				rokf=0.16;
			                          	}
						else
			                           	{
			                            	rokf=0.16*(mtk-teclmin)/(teclmax-teclmin);
			                            	}
			                	}
			                if(mpb>=p_pl_out)
			                   	{
			                       	mro*=1.0+rokf;
			                        }
			                else
			                        {
			                        mro*=(1.0+rokf*(mpb-p_ga_in)/(p_pl_out-p_ga_in));
			                        }
			                }
				    /* Coe->St transition Gerya et al., 2004, PCM */
			        p_st_in=59100.0+mtk*22.6;
			        if(mpb>p_st_in) mro*=1.06;
				    /* Pv transition, Mishin et al., 2008 with slope from Ito et al., 1990 */
			        /* Sp-out transition*/
			        p_sp_out=354000.0-mtk*40.0;
			        /* Pv-in transition*/
			        p_pv_in=352000.0-mtk*40.0;
			        if(mpb>p_pv_in)
			                {
					rokf=0.08;
			                if(mpb>=p_sp_out)
			                   	{
			                       	mro*=1.0+rokf;
			                        }
			                else
			                        {
			                        mro*=(1.0+rokf*(mpb-p_pv_in)/(p_sp_out-p_pv_in));
			                        }
			                }

			        }
				/* Ol-Sp and Pv transitions in the mantle */
				if(mm2>=9 && mm2<=14) 
			        {
				    /* Ol-Sp transition, Katsura & Ito, 1989 */
			        /* Ol-out transition*/
			        p_ol_out=91000.0+mtk*27.0;
			        /* Sp-in transition*/
			        p_sp_in=66000.0+mtk*39.0;
			        /*Limit width of Sp-Ol transition to 2 kbar */
			        if(p_sp_in>p_ol_out-2000.0) p_sp_in=p_ol_out-2000.0;
			        if(mpb>p_sp_in)
			                {
					rokf=0.06;
			                if(mpb>=p_ol_out)
			                   	{
			                       	mro*=1.0+rokf;
			                        }
			                else
			                        {
			                        mro*=(1.0+rokf*(mpb-p_sp_in)/(p_ol_out-p_sp_in));
			                        }
			                }
				    /* Pv transition, Ito et al., 1990 */
			        /* Sp-out transition*/
			        p_sp_out=304000.0-mtk*40.0;
			        /* Pv-in transition*/
			        p_pv_in=302000.0-mtk*40.0;
			        if(mpb>p_pv_in)
			                {
					rokf=0.11;
			                if(mpb>=p_sp_out)
			                   	{
			                       	mro*=1.0+rokf;
			                        }
			                else
			                        {
			                        mro*=(1.0+rokf*(mpb-p_pv_in)/(p_sp_out-p_pv_in));
			                        }
			                }
			        }
				}
				
				/* ========================== end */
			
				/* Test Heat conductivity k=ko/(1+b*(T-To)/To) */
				if (markkt[mm2]<0) mkt=-markkt[mm2]/(1.0+markkf[mm2]*(mtk-markkp[mm2])/markkp[mm2]);
				mht=markht[mm2];

				/* 4. Molten rocks */
				// Note: Calls viscalc only for melted rocks here !
				if (mm2>20) { meltpartomp(mtk,mpb,markx[mm1],marky[mm1],mm1,mm2,&mro,&mbb,&maa,&mnu,&mcp,&mkt,&mgg,&mxmelt,&mhlatent); } 

				/* ~X Thermodynamic database use for water  */
				// ATAT QUESTION TARAS: I would suggest to remove this if-statement, since it is harldy use. Do you agree? To avoid what was it included?
				/* Density && Water wt% save */
				// Hardly used; something that is larger than air or water in rock type number, but has negative density, at very start of model..	
				if ( densimod==3 && mm2>1 && (timesum<=1e+11 || markd[mm1]<=0) )
				{
					// tdbasecalc(mtk,mpb,mm2,mm1);
					// markw[mm1]=eps[42];
					// markd[mm1]=mro;
					// ATATOMP QUESTION TARAS: the above form of mro means that it comes from dencalcomp or meltpartomp, and not from tdbasecalc above that stored is as eps.. ; what you wanted ?
					// Use capital letters since Taras does not always assign eps value into this loop for m.. (eg mkt in next densimod2 call few lines below)
					// ATATOMP QUESTION TARAS: Is this intentionally? With what purpose ? eg for mkt in  next densimod2 call few lines below
					tdbasecalcomp(markx[mm1],marky[mm1],mtk,mpb,mm2,mm1,m10,&Mgg,&Mro,&Mwa,&Mcp,&Mbb,&Maa,&Mdhh,&Mkt);
					markw[mm1]=Mwa;
					markd[mm1]=mro;
				}

#elif setup<10
				// In lab have constant density per rocktype and no thermal evolution, so much faster
				mro=markro[mm2];
#endif

				/* ---> (3) Marker rheology: calculate viscosity and stresses <--- */
				mdi0=0;
				mdi1=1.0;
#if setup>9
				if(mm2<=20)
#endif 
				{
					// yn = 1 now means that plasticity IS executed in this routine call to viscalc! This is the only successfull call within current code for normal, non-melting rocks !	
					viscalcomp(mtk,mpb,markx[mm1],marky[mm1],markv[mm1],markwa[mm1],markk[mm1],markp[mm1],markt[mm1],markexx[mm1],markexy[mm1],markxx,markxy,marke,mm1,mm2,1,m10,&mnu,&mdi0);
	
					/* XXX Density correction for the dilation angle XXX = not executed */
					if(markf0[mm2]>0 && markf1[mm2]>0 && marke[mm1]>0)
					{
						/* Second invariant of viscoplastic strain calc, check */
						sigin=pow(markxx[mm1]*markxx[mm1]+markxy[mm1]*markxy[mm1],0.5);
						epsin=marke[mm1]-sigin/2.0/markgg[mm2];
						if(epsin>markf1[mm2]) epsin=markf1[mm2];
						if(epsin>0) mdi1=exp(-2.0*epsin*markf0[mm2]);
					}
				}

				msxxe=markxx[mm1];
				msxye=markxy[mm1];
				mexxe=markexx[mm1];
				mexye=markexy[mm1];
				mgg=markgg[mm2];
			
				/* Min,Max NU limitation */
				if(mnu<nubeg) mnu=nubeg; if(mnu>nuend) mnu=nuend;

				/* Water/Air account */
#if setup>9
				if(mm2<2)
				{
					markd[mm1]=mro;
					mdi0=0;
					mdi1=1.0;
				}
#endif
				/* End Water/Air account  */
			
				/* Calc log density derivative, save new density */
				if(markd[mm1]<=0 || densimod==0) 
					{markd[mm1]=mro; mdi0=0;}
				mdro=0; 
				maa=0; 
				mdi1=1.0; 
			
#if setup>9
				if(timestepe) 
				{
					mdro=mro/markd[mm1];
					mdro=log(mdro)-mdi0;
					//if(epsin>0 && debugmod) {fprintf(fp_log,"d %ld %d  %e %e   %e %e   %e %e   %e %e %e %e",mm1,mm2,markx[mm1],marky[mm1],marke[mm1]*2.0*markgg[mm2],sigin,marke[mm1],epsin,-2.0*epsin*markf0[mm2],markd[mm1],mro,mdro);getchar();}
				}
#endif
			
				/* Save new density */
				mdro=-mdi0; 
				markd[mm1]=mro;
				/* Correct new density for dilation */
				mro*=mdi1;
			
				/* Saving marker viscosity */
				markv[mm1]=mnu;
			
				// if(debugmod) {fprintf(fp_log,"num=%ld type=%d  x=%e y=%e mpb=%e mtk=%e nu=%e ro=%e cp=%e kt=%e ht=%e",mm1,mm2,markx[mm1],marky[mm1],mpb,mtk,mnu,mro,mcp,mkt,mht);getchar()};
			
				/* --> (4) Interpolation from markers to 4 corners of the cell ====================================*/
				
				/* Marker weight calculation using dimension of current Cell */
				celdx=gx[m10+1]-gx[m10];
				celdy=gy[m20+1]-gy[m20];
				swt1=1.0/celdx/celdy;
				/* Marker weights calculation using dimension of current Cell */
				celdx=(markx[mm1]-gx[m10])/(gx[m10+1]-gx[m10]);
				celdy=(marky[mm1]-gy[m20])/(gy[m20+1]-gy[m20]);
				if (celdx<0 || celdy<0 || celdx>1.0 ||celdy>1.0) {fprintf(fp_log," WARNING !!! num=%ld type=%d  x=%e y=%e celdx=%e celdy=%e",mm1,mm2,markx[mm1],marky[mm1],celdx,celdy); fflush(fp_log); getchar();}
			
				/* --- Interpolate ro,nu etc to nodes using interpolation coefficients --- */
				for (m1=0;m1<4;m1++)
				{
					/* Marker weight calculation using dimension of current Cell */
					/* Different corners */
					/* 0  2 */
					/* 1  3 */
					switch(m1)
					{
						case  0: 
						/* Calc node number */
						m3=m10*ynumy+m20;
						
						/* Add shear viscosity Nu */
						if (celdx<0.5 && celdy<0.5) 
						{
							dx=1.0-2.0*celdx;
							dy=1.0-2.0*celdy;
							swt=swt1*dx*dy;
							Nu0[m3]+=mnu*swt;
							Gg0[m3]+=mnu/mgg*swt;
							Sxye0[m3]+=msxye*swt;
							Exye0[m3]+=mexye*swt;
							Sol0[nodenum+m3]+=swt;
						}
						/* Add Vx and Mx from markers */
						if (celdx<0.5) 
						{
							dx=1.0-celdx;
							dy=1.0-ABSV(celdy-0.5);
							swt=swt1*dx*dy;
							Mvx0[m3]+=markvx[mm1]*mro*swt;
							Mrx0[m3]+=mro*swt;
							Sol0[nodenum2+m3]+=swt;
						}
						/* Add Vy and My from markers */
						if (celdy<0.5) 
						{
							dx=1.0-ABSV(celdx-0.5);
							dy=1.0-celdy;
							swt=swt1*dx*dy;
							Mvy0[m3]+=markvy[mm1]*mro*swt;
							Mry0[m3]+=mro*swt;
							Sol1[nodenum2+m3]+=swt;
						}
						
						// Calculate standard weight for physical properties
						swt=swt1*(1.0-celdx)*(1.0-celdy); 
						break;
						
						case  1: 
						/* Calc node number */
						m3=m10*ynumy+m20+1;
						/* Add shear viscosity Nu */
						if (celdx<0.5 && celdy>0.5) 
						{
							dx=1.0-2.0*celdx;
							dy=2.0*celdy-1.0;
							swt=swt1*dx*dy;
							Nu0[m3]+=mnu*swt;
							Gg0[m3]+=mnu/mgg*swt;
							Sxye0[m3]+=msxye*swt;
							Exye0[m3]+=mexye*swt;
							Sol0[nodenum+m3]+=swt;
						}
						/* Add Vy and My from markers */
						if (celdy>0.5) 
						{
							dx=1.0-ABSV(celdx-0.5);
							dy=celdy;
							swt=swt1*dx*dy;
							Mvy0[m3]+=markvy[mm1]*mro*swt;
							Mry0[m3]+=mro*swt;
							Sol1[nodenum2+m3]+=swt;
						}
						
						// Calculate standard weight for physical properties
						swt=swt1*(1.0-celdx)*celdy; 
						break;
						
						case  2: 
						/* Calc node number */
						m3=(m10+1)*ynumy+m20;
						/* Add shear viscosity Nu, Sxy */
						if (celdx>0.5 && celdy<0.5) 
						{
							dx=2.0*celdx-1.0;
							dy=1.0-2.0*celdy;
							swt=swt1*dx*dy;
							Nu0[m3]+=mnu*swt;
							Gg0[m3]+=mnu/mgg*swt;
							Sxye0[m3]+=msxye*swt;
							Exye0[m3]+=mexye*swt;
							Sol0[nodenum+m3]+=swt;
						}
						/* Add Vx and Mx from markers */
						if (celdx>0.5) 
						{
							dx=celdx;
							dy=1.0-ABSV(celdy-0.5);
							swt=swt1*dx*dy;
							Mvx0[m3]+=markvx[mm1]*mro*swt;
							Mrx0[m3]+=mro*swt;
							Sol0[nodenum2+m3]+=swt;
						}
						
						// Calculate standard weight for physical properties
						swt=swt1*celdx*(1.0-celdy); 
						break;
						
						case  3: 
						/* Calc node number */
						m3=(m10+1)*ynumy+m20+1;
						/* Add shear viscosity Nu */
						if (celdx>0.5 && celdy>0.5) 
						{
							dx=2.0*celdx-1.0;
							dy=2.0*celdy-1.0;
							swt=swt1*dx*dy;
							Nu0[m3]+=mnu*swt;
							Gg0[m3]+=mnu/mgg*swt;
							Sxye0[m3]+=msxye*swt;
							Exye0[m3]+=mexye*swt;
							Sol0[nodenum+m3]+=swt;
						}
						
						// Add values to central node once
						// Brackets determine the scope of the to here limited variables, but currently make no difference
						{
							dx=1.0-2.0*ABSV(celdx-0.5);
							dy=1.0-2.0*ABSV(celdy-0.5);
							swt=swt1*dx*dy;
							Nd0[m3]+=mnu*swt;
							Gd0[m3]+=mnu/mgg*swt;
							Sxxe0[m3]+=msxxe*swt;
							Sppe0[m3]+=markp[mm1]*swt;
							Sbritn0[m3]+=msbrit[mm1]*swt; // Yield stress in the pressure node
							Exxe0[m3]+=mexxe*swt;
							Dro0[m3]+=mdro*swt;
							Drp0[m3]+=maa*swt;
							Sol1[nodenum+m3]+=swt;
						}
						
						// Calculate standard weight for physical properties
						swt=swt1*celdx*celdy; 
						break;
					}
					// End switch of weight calculation
				
					/* Add Physical Properties: ro,nu, etc. */
					// fprintf(fp_log,"num=%ld type=%d x=%e y=%e cell=%ld dx=%e dy=%e swt=%e",mm1,mm2,markx[mm1],marky[mm1],m3,dx,dy,swt);getchar();
					// nu0[m3]+=mnu*swt;
				
					// ATAT TARAS	 Why use mcp*MRO? in routines calculate purely from markcp... not done in manuele's, later in heat mrocp, but unrelated
					Ro0[m3]+=mro*swt;
					Et0[m3]+=mbb*swt;
					Cp0[m3]+=mcp*mro*swt;
					Kt0[m3]+=mkt*swt;
					Ht0[m3]+=mht*swt;
					Sol0[m3]+=swt;
				
					/* Add T */
					if(!markim[mm2]) 
					{
						Tk[m3]+=mtk*swt;
						Sol1[m3]+=swt;
					}
				}
				/* End Interpolation from markers to nodes ====================================*/
			}
		}
		
		// Add interpolation arrays from different processors and free their memory
#pragma omp critical (sumsolarrays)
		{
			for (m3=0;m3<nodenum;m3++)
			{
				nu0[m3]+=Nu0[m3];
				nd0[m3]+=Nd0[m3];
				gd0[m3]+=Gd0[m3];
				gg0[m3]+=Gg0[m3];
				ro0[m3]+=Ro0[m3];
				cp0[m3]+=Cp0[m3];
				kt0[m3]+=Kt0[m3];
				ht0[m3]+=Ht0[m3];
				dro0[m3]+=Dro0[m3];
				drp0[m3]+=Drp0[m3];
				sxxe0[m3]+=Sxxe0[m3];
				sxye0[m3]+=Sxye0[m3];
				sppe0[m3]+=Sppe0[m3];
				sbritn0[m3]+=Sbritn0[m3];
				exxe0[m3]+=Exxe0[m3];
				exye0[m3]+=Exye0[m3];
				et0[m3]+=Et0[m3];
				tk0[m3]+=Tk[m3];
				mrx0[m3]+=Mrx0[m3];
				mry0[m3]+=Mry0[m3];
				mvx0[m3]+=Mvx0[m3];
				mvy0[m3]+=Mvy0[m3];
				sol0[m3]+=Sol0[m3];
				sol1[m3]+=Sol1[m3];
				sol0[nodenum+m3]+=Sol0[nodenum+m3];
				sol1[nodenum+m3]+=Sol1[nodenum+m3];
				sol0[nodenum2+m3]+=Sol0[nodenum2+m3];
				sol1[nodenum2+m3]+=Sol1[nodenum2+m3];		   
			}	
		}
						
		// Free dynamically allocated interpolation arrays	
		free(Nu0);
		free(Nd0); 
		free(Gg0); 	
		free(Gd0);
		free(Ro0); 
		free(Sxxe0); 
		free(Sppe0); 
		free(Sbritn0);	
		free(Sxye0); 
		free(Exxe0); 
		free(Exye0); 
		free(Et0);
		free(Dro0); 	
		free(Drp0); 
		free(Cp0); 
		free(Kt0); 
		free(Ht0); 	
		free(Tk); 	
		free(Mrx0); 
		free(Mry0); 
		free(Mvx0); 	
		free(Mvy0); 
		free(Sol0);
		free(Sol1); 
	}
	// End OMP section marker to node interpolation

#if setup>9
	if (n0==1){fclose(flfric);}
#endif

	if (printmod==10000) fprintf(fp_log,"\n Time taken for rho and vis calc + M->N1 in ronurecalc = %e s \n",omp_get_wtime()-start);
	start=omp_get_wtime();
	
	/* Recalculate ro[] nu[] */
	for (m1=0;m1<xnumx;m1++)
	{
		for (m2=0;m2<ynumy;m2++)
		{
			/* Current node num, wt */
			m3=m1*ynumy+m2;
			
			/* Shear viscosity recalc check */
			if(sol0[nodenum+m3])
			{
				// Boundary Condition Viscosity (set in mu)
				if(mu[m3] && (timesum<timebond || m1<=2 || m2<=2 || m1>=xnumx-4 || m2>=ynumy-3)) 
				{
					// BC value defined in init.t3c
					if(mu[m3]>0)
					{
						nu0[m3]=mu[m3];
					}
					else
					{
						nu0[m3]/=sol0[nodenum+m3];
						if(nu0[m3]>-mu[m3]) nu0[m3]=-mu[m3];
					}
				} 
				// Rest; solution
				else
				{
					nu0[m3]/=sol0[nodenum+m3];
				} 
				/* Min,Max NU limitation */
				if(nu0[m3]<nubeg) nu0[m3]=nubeg; if(nu0[m3]>nuend) nu0[m3]=nuend;
				/* Min,Max NU definition for nu contrast limit */
				if(numin==0 || nu0[m3]<numin) numin=nu0[m3]; if(numax==0 || nu0[m3]>numax) numax=nu0[m3];
				nu[m3]=nu0[m3];
					
				/* Elastic shear stress Sxy recalc */
				sxye[m3]=sxye0[m3]/sol0[nodenum+m3];
				exye[m3]=exye0[m3]/sol0[nodenum+m3];
					
				/* Shear shear modulus recalc */
				gg[m3]=nu[m3]/(gg0[m3]/sol0[nodenum+m3]);
					
				/* Reset weight */
				sol0[nodenum+m3]=0;
			} 
				
			/* Normal viscosity recalc check */
			if(sol1[nodenum+m3])
			{
				if(mu[m3] && (timesum<timebond || m1<=2 || m2<=2 || m1>=xnumx-4 || m2>=ynumy-3)) 
				{
					if(mu[m3]>0)
					{
						nd0[m3]=mu[m3];
					}
					else
					{
						nd0[m3]/=sol1[nodenum+m3];
						if(nd0[m3]>-mu[m3]) nd0[m3]=-mu[m3];
					}
				} 
				else
				{
					nd0[m3]/=sol1[nodenum+m3];
				} 
				/* Min,Max NU limitation */
				if(nd0[m3]<nubeg) nd0[m3]=nubeg; if(nd0[m3]>nuend) nd0[m3]=nuend;
				/* Min,Max NU definition for nu contrast limit */
				if(numin==0 || nd0[m3]<numin) numin=nd0[m3]; if(numax==0 || nd0[m3]>numax) numax=nd0[m3];
				nd[m3]=nd0[m3];
					
				/* Elastic Normal stress recalc */
				sxxe[m3]=sxxe0[m3]/sol1[nodenum+m3];
				sppe[m3]=sppe0[m3]/sol1[nodenum+m3];
				sbritn[m3]=sbritn0[m3]/sol1[nodenum+m3];
				exxe[m3]=exxe0[m3]/sol1[nodenum+m3];
				/* Density changes recalc */
				dro[m3]=dro0[m3]/sol1[nodenum+m3];
				drp[m3]=drp0[m3]/sol1[nodenum+m3];
					
				/* Normal shear modulus recalc */
				gd[m3]=nd[m3]/(gd0[m3]/sol1[nodenum+m3]);
					
				/* Reset weight */
				sol1[nodenum+m3]=0;
			} 
				
			/*  Vx Mx recalc check */
			if(sol0[nodenum2+m3])
			{
				/* Material constants recalc */
				mvx[m3]=mvx0[m3]/mrx0[m3];
				mrx[m3]=mrx0[m3]/sol0[nodenum2+m3];
				sol0[nodenum2+m3]=0;
			} 
				
			/*  Vy My recalc check */
			if(sol1[nodenum2+m3])
			{
				/* Material constants recalc */
				mvy[m3]=mvy0[m3]/mry0[m3];
				mry[m3]=mry0[m3]/sol1[nodenum2+m3];
				sol1[nodenum2+m3]=0;
			} 
				
			/* Other variables recalc check */
			if(sol0[m3])
			{
				/* Material constants recalc */
				ro[m3]=ro0[m3]/sol0[m3];
#if setup>9
				if(gy[m2]<waterlev && ro[m3]<1000.1) ro[m3]=1.0;
				if(gy[m2]>=waterlev && ro[m3]<1000.1) ro[m3]=1000.0;
#endif
				et[m3]=et0[m3]/sol0[m3];
				cp[m3]=(cp0[m3]/sol0[m3])/ro[m3];
				kt[m3]=kt0[m3]/sol0[m3];
				ht[m3]=ht0[m3]/sol0[m3];
					
				/* Advective addition for T K in nodes recalc */
				if (sol1[m3]) 
				{
					tk[m3]=tk0[m3]/sol1[m3]; 
					sol1[m3]=0;
				}
					
				/* Reset weight */
				sol0[m3]=0;
			}
		}
	}
		
	if (printmod) fprintf(fp_log,"Min, Max viscosity %e %e \n",numin,numax); fflush(fp_log);

	/* Reset advective temperature */
	for (m3=0;m3<nodenum;m3++) {tk3[m3]=0;} 
	
	/* Check Upper/Lower limits for nu[] after given contrast */
	if(nucontr>1.0 && numin>0) numax=numin*nucontr;
	if(nucontr<1.0 && numax>0) numin=numax*nucontr;
	for (m3=0;m3<nodenum;m3++)
	{
		if(nu[m3]<numin) nu[m3]=numin; if(nu[m3]>numax) nu[m3]=numax;
		if(nd[m3]<numin) nd[m3]=numin; if(nd[m3]>numax) nd[m3]=numax;
	}
	
	/* Water/air density */
#if setup>9
	for (m1=0;m1<xnumx;m1++)
		for (m2=0;m2<ynumy;m2++)
	{
		m3=m1*ynumy+m2;
		if(gy[m2]<waterlev && ro[m3]<1000.1) ro[m3]=1.0;
		if(gy[m2]>=waterlev && ro[m3]<1000.1) ro[m3]=1000.0;
	}
#endif
	
	/* ---> 5. Set Boundary conditions for T <---*/
	if (printmod) fprintf(fp_log,"\n AVERAGE TEMPERATURE CORRECTION FOR BOUNDARY CONDITIONS ...\n"); fflush(fp_log);
	tkrecalc();
	if (printmod) fprintf(fp_log,"AVERAGE TEMPERATURE OK!\n"); fflush(fp_log);
	/* Adiabate computing */
	if(1==0 && timesum<3.15576e+7*1e+3) 
	{
		/* Lower boundary TK - Node Cycle */
		for (m1=0;m1<xnumx;m1++)
		{
			/* Cur Line Num in bondm[] */
			m2=(m1+1)*ynumy-1;
			m3=bondm[m2+nodenum3];
			if(m3) 
				{bondv[m3][0]=tk[m2-1]*2.0-tk[m2-2];}
		}
	}
	
	if (printmod==10000) fprintf(fp_log,"\n Time taken for M->N2 in ronurecalc = %e s \n",omp_get_wtime()-start);
	
	/* END ADD MARKERS TO THE v-CELLS ========================== */
}
/* End ro[],nu[] recalc after marker positions - routine */


/* Calc density for given P,T */
void dencalcomp(double mtk, double mpb, double x, double y, int mm2, double *mro, double *mbb, double *maa)
	/* mtk - T, K */
	/* mpb - P, bar */
	/* x,y - XY location of point for Vx,Vy calc */
	/* mm2 - Rock number */
{
	/* Ro=ro0*(1-bro*(TK-298.15))*(1+aro*(Pkbar-0.001)) */
	*mro=markro[mm2]*(1.0-markbb[mm2]*(mtk-298.15))*(1.0+markaa[mm2]*(mpb-1.0)*1e-3);
	/* Adiabatic term: al=bro/(1-bro*(Tk-298.15)) */
	*mbb=markbb[mm2]/(1.0-markbb[mm2]*(mtk-298.15));
	/* Compressibility: be=aro/(1+aro*(Pkbar-0.0001) */
	*maa=1.e-8*markaa[mm2]/(1.0+markaa[mm2]*(mpb-1.0)*1e-3);
	
	/* Constant density */
	if (densimod==0) *mro=markro[mm2];
}
/* End OMP Calc density for given P,T */


/* OMP Antigorite weakening of mantle */
void antigoromp(double mtk, double mpb, double x, double y, long int mm1, long int m10, char markt[])
	/* mtk - T, K */
	/* mpb - P, bar */
	/* x,y - XY location of point for Vx,Vy calc */
	/* mm1 - mark number */
	/* m10 - Up Left Node X,Y Num */
{
	/* Val buffer */
	double k1,sy1,e,hydry,yfiltr,hydryl,tsubd,vxs,vys;

	/* Check marker type */
	if(markt[mm1]!=11 && markt[mm1]!=13) return;

	/* Relativ Normalized coord Calc */
	e=(x-gx[m10])/(gx[m10+1]-gx[m10]);
	/* Erosion surface; oceanic crust top */
	sy1=(e*ep[m10+1]+(1.0-e)*ep[m10]);

	/* Antigorite weakening of mantle above oceanic crust */
	/* Atg stability field after Schmidt and Poli, 1998 */
	if((y-sy1)>63000.0)
	{
		k1=1013.17699-0.060387633e-3*(y-sy1)-0.004289442e-6*(y-sy1)*(y-sy1);
	}
	else
	{
		k1=751.490422+6.00773668e-3*(y-sy1)-0.034690759e-6*(y-sy1)*(y-sy1);
	}
	
	/* Change marker Type */
	/* Serpentinized (13) - to hydrated (11) */
	if(k1<=mtk && markt[mm1]==13) markt[mm1]=11;
	/* Hydrated(11) -  to serpentinized (13) */
	if(k1>mtk && markt[mm1]==11) markt[mm1]=13;
}
/* OMP End Antigorite weakening of mantle */



/* Nu calc after reological equation */
// Uses timestepe or computational visco-elastic timestep
/* P-T-stress dependent rheology without/with brittle/ductile transition */
/* Reological equations */
/* Stress>SScr */
/* Power law dislocation creep: SSii={NU0*EEii*exp[(E+PV)/RT]}^(1/n) */
/* Effective viscosity: NU=1/2*{NU0*exp[(E+PV)/RT]}^(1/n)*EEii^[(1-n)/n] */
/* Stress<SScr */
/* Newtonian diffusion   creep: SSii=NU1*EEii*exp[(E+PV)/RT] */
/* Effective viscosity: NU=NU0/2*exp[(E+PV)/RT] */
/* NU1=NU0/SScr^(n-1) */
/* SScr - dislocation, diffusion transition stress */
/* SSii - second invariant of deviatoric stress tensor */
/* EEii - epsin - second invariant of strain rate tensor */
/* E - activation energy, J */
/* V - activation volume, J/bar */
/* R - gase constant 8.314 J/K */
/* Viscosity NU  calc after reological equations */
/* NU=SSii/(2*EEii) */
/* Brittle - Ductile transition */
/* sbrit=MINV(0.85e+5*pb,60e+6+0.6e+5*pb)*lambda;  (Schott & Schmeling, 1998) */
/* sbrit=MINV(0.667e+5*pb,51.2e+6+0.512e+5*pb)*lambda; (Brace & Kohlsstedt, 1980) */
void viscalcomp(double mtk, double mpb, double cmx, double cmy, double Markv, double Markwa, double Markk, double Markp, double Markt, double Markexx, double Markexy, double Markxx[], double Markxy[], double Marke[], long int mm1, int mm2, int yn, long int m10,double *mnu, double *mdi0)
	/* mtk - T, K */
	/* mpb - P, bar */
	/* cmx,cmy - XY location of point for Vx,Vy calc */
	/* mm1 - Marker number */
	/* mm2 - rock type */
	/* yn - plastic reset yes(1)/no(0) - switch from version 1 to 2 ! */
	// bbrit_cor - slip velocity dependent correction for friction coefficient 
{
	/* Val buffer */
	double xnu,nnu,e,n,rt=8.314*mtk,k1,e1,epsin,sigin,sduct,sbrit,nueff,strain,abrit,bbrit,nubrit,nunewt,nupowl,nuduct;
	/* Reological Eq par */
	double sy1,lamb,xelvis,sxxnew,sxynew,siginnew,mnu0,mnu1,mnu2,siginnew0,siginnew1,siginnew2,dsiginnew0,dsiginnew1,dsiginnew2;
	/* Counters */
	long int m1;
	int ncount=0; 
	// Slip velocity dependent friction
	double relvw=60,dvw,dvs;  
	
	/* Melted rocks */
	// But incoming mm2 is already mm2_actual-20
#if setup>9
	if (mm2>20) { *mnu = markn0[mm2]; return; }
#endif
	
	/* Non-melted rocks, mm2 <= 20 */
	/* Calc effective strain rate, stress after second strain rate Tenzor invariant EEii=(1/2SUM(EPSik^2))^(1/2) */
	// Interpolation to these markers done at end of viterate() of previous timestep
	epsin=pow(Markexx*Markexx+Markexy*Markexy,0.5);
	sigin=pow(Markxx[mm1]*Markxx[mm1]+Markxy[mm1]*Markxy[mm1],0.5);


	/* --- 1. Calculate components of brittle strength; cohesion, friction, and hydrostatic pore pressure weakening factor --- */

	// - Lambda brittle weakening factor for hydrostatic pore pressure -

	/* Up Left Node X,Y Num */
	m1=m10;
	
	/* Relative normalized coord calc */
	e=(cmx-gx[m1])/(gx[m1+1]-gx[m1]);
	n=(e*ep[m1+1]+(1.0-e)*ep[m1]);
	
	// Pore fluid pressure correction: lamb = Pf/Ps
	lamb=markll[mm2]; 
#if setup>9
	// Predefine fluid pressures near surface 
	if ((cmy-n)<=0) lamb=hidrl;
	if ((cmy-n)>0 && (cmy-n)<hidry) lamb=hidrl*(1.0-(cmy-n)/hidry)+lamb*(cmy-n)/hidry;

	// Lower friction in fluid/melt present areas
	if (Markwa==1) {lamb=lambfld;}
#endif
	
	/* - Strain weakening - */

	strain=Marke[mm1];

	/* A,B coefficients calc depending on integral strain */

        abrit=marka0[mm2];
        bbrit=markb0[mm2];
        if(strain>marke1[mm2])
        {
                abrit=marka1[mm2];
                bbrit=markb1[mm2];
        }
        else
        {
                if(strain>marke0[mm2] && marke1[mm2]>marke0[mm2])
                {
                        abrit=marka0[mm2]+(marka1[mm2]-marka0[mm2])*(strain-marke0[mm2])/(marke1[mm2]-marke0[mm2]);
                        bbrit=markb0[mm2]+(markb1[mm2]-markb0[mm2])*(strain-marke0[mm2])/(marke1[mm2]-marke0[mm2]);
                }
        }
	
	/* --- End calculation of brittle strength components; cohesion, friction, and hydrostatic pore pressure weakening factor --- */


	/* --- Start ductile viscosity calculation -------------------------------------------*/
	/* Inverted value of newtonian NU set */
	nunewt=0;

	/* Inverted value of power-low NU set */
	nupowl=0;

	/* Check for the presence of ductile rheology */
	// For more viscosity options, see codes of version 1
	if (marknu[mm2])
	{
		/* A)  Simple Newtonian rheology */
		// - used in laboratory model of van Dinther et al., JGR, 2013a -
		/* Newtonian creep: SSii=NU0*2.0*EEii */
		/* Effective viscosity: NU=NU0 */
		/* Effective viscosity member in Stoks: NUs=NU */
		if(markdh[mm2]==0 && markdv[mm2]==0 && (markss[mm2]==0 || markmm[mm2]==1.0))
		{
			/* Inverted value of newtonian NU calc */
			nunewt=1.0/marknu[mm2];
		}

		/* --> D)  P-T-stress dependent rheology without/with brittle/ductile transition <--*/
		// - used in large-scale models PhD thesis van Dinther -
		/* Reological equations */
		/* Stress>SScr */
		/* Power law dislocation creep: SSii={NU0*EEii*exp[(E+PV)/RT]}^(1/n) */
		/* Effective viscosity: NU=1/2*{NU0*exp[(E+PV)/RT]}^(1/n)*EEii^[(1-n)/n] */
		/* Effective viscosity member in Stoks: NUs=NU/n */
		/* Stress<SScr */
		/* Newtonian diffusion   creep: SSii=NU1*EEii*exp[(E+PV)/RT] */
		/* Effective viscosity: NU=NU0/2*exp[(E+PV)/RT] */
		/* Effective viscosity member in Stoks: NUs=NU */
		/* NU1=NU0/SScr^(n-1) */
		if(marknu[mm2]>0 && (markdh[mm2]!=0 || markdv[mm2]!=0) && markss[mm2]!=0 && markmm[mm2]!=1.0)
		{
			// ---> 2. Calculate ductile viscosity <--- 
			/* T-P exponent for effective NU calc */
			e1=(markdh[mm2]+markdv[mm2]*mpb)/rt;
			if(e1>150.0) e1=150.0;
			e1=exp(e1);

			/* Koef for stress independent creep NU1 calc */
			k1=marknu[mm2]/pow(markss[mm2],markmm[mm2]-1.0);

			/* Inverted value of newtonian NU calc for diffusion creep */
			nunewt=1.0/(0.5*k1*e1);
			mnu2=nunewt;

			/* Effective viscosity1 calc */
			siginnew1=siginnew=sigin;
			nupowl=0;
			
			
			// Calculate dislocation creep viscosity
			if (siginnew>0) nupowl=1.0/(0.5*siginnew*marknu[mm2]*e1/pow(siginnew,markmm[mm2]));
			mnu1=nupowl;
			//Take arithmetic average of dislocation and diffusionc creep for effective ductile viscosity
			mnu0=1.0/(mnu1+mnu2);
			
			// ---> 3. Include elastic part for estimation future viscoelastic stresses <---
			// Calculate visco-elasticity factor
			xelvis=markgg[mm2]*timestepe/(markgg[mm2]*timestepe+mnu0);
			// Calculate viscoelastic stress
			siginnew2=2.0*mnu0*epsin*xelvis+sigin*(1.0-xelvis);
			dsiginnew1=siginnew2-siginnew1;
			
			
			/* Effective viscosity2 calc */
			// See above for description. Repeated here 
			siginnew=siginnew2;
			nupowl=0;
			if (siginnew>0) nupowl=1.0/(0.5*siginnew*marknu[mm2]*e1/pow(siginnew,markmm[mm2]));
			mnu1=nupowl;
			mnu0=1.0/(mnu1+mnu2);
			// Calculate visco-elasticity factor
			xelvis=markgg[mm2]*timestepe/(markgg[mm2]*timestepe+mnu0);
			// Calculate viscoelastic stress
			siginnew=2.0*mnu0*epsin*xelvis+sigin*(1.0-xelvis);
			dsiginnew2=siginnew-siginnew2;
			
			
			/* ---> 4. Local iterations for dislocation viscosity calculation by Bisection method <--- */
			ncount=0;
			// Locally iterate over nupowl and siginnew until siginnew-siginnew0<10 Pa (or 100 iterations)
			do
			{
				// Check to prevent num issue when stress is not changing: only not true and in if sigma_2_1 = sigma_1_1 = sigma_0 : almost never no stress change?
				dsiginnew0=ABSV(dsiginnew1)+ABSV(dsiginnew2);
				if(dsiginnew0>0)
				{
					// Weigth factor: 0.5 = midpoint
					dsiginnew0=0.5;
					
					// Calculate midpoint = new estimate stress
					siginnew0=siginnew=siginnew1*(1.0-dsiginnew0)+siginnew2*dsiginnew0;

					// Update viscosity with that new stress
					nupowl=0;
					if (siginnew>0) nupowl=1.0/(0.5*siginnew*marknu[mm2]*e1/pow(siginnew,markmm[mm2]));
					mnu1=nupowl;
					mnu0=1.0/(mnu1+mnu2);
					
					// Update stress estimate with new viscosity
					xelvis=markgg[mm2]*timestepe/(markgg[mm2]*timestepe+mnu0);
					siginnew=2.0*mnu0*epsin*xelvis+sigin*(1.0-xelvis);
					
					// Calculate difference new and last stress estimate -> converging?
					dsiginnew0=siginnew-siginnew0;
					
					// Use this newest estimate for stress change to see if in same direction
					// If yes; keep going in that direction; leave oldest estimate behind 
					if((dsiginnew0>=0 && dsiginnew1>=0) || (dsiginnew0<0 && dsiginnew1<0)) 
					{
						siginnew1=siginnew0;
						dsiginnew1=dsiginnew0;
					}
					// If in opposite direction; passed optimal stress so turn back; leave newest 1 estimate behind
					else
					{
						siginnew2=siginnew0;
						dsiginnew2=dsiginnew0;
					}
				}
				ncount++;
			}
			while(ABSV(dsiginnew0)>10.0 && ncount<101);
		}
	}
	/* --- End Ductile viscosity calculation -------------------------------------------*/

	// Check ductile effective viscosity calculation
	nueff=1.0/(nunewt+nupowl);
	
	/* Mantle viscosity */	
#if setup > 9
	if((Markt==9 || Markt==10) && timesum<3.15576e+7*1e+4 && nueff<1e+20) nueff=1e+20;
#endif
	
	if(nueff<nubeg) nueff=nubeg; if(nueff>nuend) nueff=nuend;
	if(nueff<markn0[mm2]) nueff=markn0[mm2]; if(nueff>markn1[mm2]) nueff=markn1[mm2];
	nuduct=nueff;
	*mdi0=0;
	
	
	
	/* ------------------ Calculate viscoplastic viscosity ---------------------------- */
	
	// Calculate brittle strength - sbrit -
	// Plasticity switched off when both terms in the yield strength formulation are 0
	if(((1-markll[mm2])*markb0[mm2] || abrit) && epsin)
	{		
		// --- Strong slip velocity dependency of friction coefficient --- 
		// After Burridge and Knopoff (1967), Ampuero and Ben-Zion (2008), etc.
		// Adapt friction parameters based on x-location (lab) or temperature (large-scale)
#if setup==10
		if (mm2==7 && cmx>(700e3-shift_km) && cmx<(1150e3-shift_km) && cmy<100e3 && veldepfric==1)
#endif

#if setup==11	// Including off-megathrust rate weakening
		if (cmx>(700e3-shift_km) && cmx<(1150e3-shift_km) && cmy<100e3 && veldepfric==1)
#endif

#if setup==12   // Collisional setup - L. Dal Zilio
                if (cmx>(1700e3-shift_km) && cmx<(2150e3-shift_km) && cmy<100e3 && veldepfric==1)
#endif

#if setup < 10
		if (mm2==5 && veldepfric==1)
#endif
		{
			// Calculate Relative amount of Velocity-Weakening vs Velocity-Strengthening
#if setup>9
			// Velocity-strengthening region 
			if (Markk<=tk_updipsez0) // && mm2==7)
				{ relvw = 0; }
			// Transitions to seismogenic zone: updip
			else if (Markk>tk_updipsez0 && Markk<tk_updipsez1) // && mm2==7)
				{ relvw = (Markk-tk_updipsez0)/(tk_updipsez1-tk_updipsez0);	}
			// Velocity-weakening for Seismogenic Zone (and off-megathrust region)
			else
				{ relvw = 1;}
			// Note for the off-events setup there is no strengthening outside the subduction channel of basaltic crust
	
#elif setup < 10
			// Change mm2 locally in this viscalc-routine, so that also for viscosity, shear modulus, Pf/Ps etc use this 
			// Seismogenic Zone = velocity-weakening
			if (cmx >= end_updip && cmx <= start_downdip)
				{ relvw = 1; mm2 = 6; }
			// Transitions to seismogenic zone: updip
			else if (cmx >= start_updip && cmx <= end_updip)
				{ relvw = (cmx-start_updip)/(2*half_range);}
			// Transitions away from seismogenic zone: downdip
			else if (cmx >= start_downdip && cmx <= end_downdip)
				{ relvw = 1 - ( cmx-start_downdip )/( 2*half_range );}
			// Velocity-strenghtening region
			else
				{ relvw = 0; }
#endif			

			// Calculate slip-rate dependent change of coefficients
			mvslip[mm1] = 2.0*epsin*res_high;
			dvw = (1-mgamma_vw)+mgamma_vw/(1.0+mvslip[mm1]/mvc_vw);
			dvs = (1-mgamma_vs)+mgamma_vs/(1.0+mvslip[mm1]/mvc_vs);
       	
			// Change friction coefficient accordingly
			bbrit = mus_vs*dvs + relvw*(markb0[mm2]*dvw-mus_vs*dvs);
			
			// Change cohesion as a function of slip velocity, if desired (if marka0[5]=~marka0[6])
#if setup>9
			abrit = marka0[mm2];
#elif setup < 10
			abrit = marka0[5] + relvw*(marka0[6]-marka0[5]);
#endif		
        
			// Iterate locally to obtain stable estimate of slip-rate
			sbrit=abrit+bbrit*(1-lamb)*Markp;
			if(sbrit>0 && Markv>0)
			{
				for(ncount=0;ncount<5;ncount++)
				{
					if(sbrit>0)
					{
						// epsin = 0.5* sbrit/eta in viscous formulation
						mvslip[mm1] = sbrit/Markv*res_high;
						dvw = (1-mgamma_vw)+mgamma_vw/(1.0+mvslip[mm1]/mvc_vw);
						dvs = (1-mgamma_vs)+mgamma_vs/(1.0+mvslip[mm1]/mvc_vs);
                  
						bbrit = mus_vs*dvs + relvw*(markb0[mm2]*dvw-mus_vs*dvs);
                  
						sbrit=abrit+bbrit*(1-lamb)*Markp;
					}
					else
					{
						fprintf(fp_log,"LOOK: Sbrit is <= 0 within v-w loop: %e, abrit = %e, bbrit = %e, pr = %e, markvis = %e, x = %e, y = %e \n",sbrit,abrit,bbrit,Markp,Markv,cmx,cmy); fflush(fp_log);
					}
				}
			}
			
			// Calculate average value stresses and strainrates seismogenic zone
			// But here do not have proper stress and vel(e) yet ! Only yield strength ..
#if setup < 10
			if (relvw==1 && cmy<=gy[n_glayer+1])
			{
				sbrit_ave = sbrit_ave + (abrit + bbrit*(1-lamb)*Markp);
				count_sezm = count_sezm + 1;	
			}
#endif  
		}
		// In case of no rate dependency also calculate yield strength
		else
		{
			sbrit=abrit+bbrit*(1-lamb)*Markp;
		}

		// Check strength values
		if(sbrit<0) sbrit=0;
		if(sbrit>marks1[mm2]) sbrit=marks1[mm2];
		
		// Save frictional properties to file for analyses
		// Save time and space by only doing at last timestep in prn output cycle
#if setup > 9
		if (n0==1 && start_cond==1 && relvw<=1.0)
			{ fprintf(flfric," %d %e %e %e %e %e %e %e %e \n", mm2, cmx, cmy, Markp, sbrit, Markxx[mm1], Markxy[mm1],lamb,bbrit); }
#endif


		// Store yield stress for post-processing and interpolation to nodes
		msbrit[mm1]   = sbrit;
		// Store old-stress
		msii_old[mm1] = sigin;
		
		/* ---> 5. ! Viscoelastic case ! <--- */
		if(stoksmod && timestepe && epsin)
		{
			/* Future plastic creep */
			/* Future stresses calc */
			xelvis=markgg[mm2]*timestepe/(markgg[mm2]*timestepe+nueff);
			siginnew=2.0*nueff*epsin*xelvis+sigin*(1.0-xelvis);
			
			// Plastic yielding if new estimate or stress of previous timestep exceeds strength	
			if(sbrit<siginnew || sbrit<sigin)
			{					
				/* Executing plasticity by reseting stresses and viscosities */
				// Note yn is defined at call to this viscalc routine 
				if(yn==1)
				{
					/* XXX Density correction for the dilation angle XXX */
					// We do not use dilation !	Not for any rock type
					if(markf0[mm2]>0 && markf1[mm2]>0)
					{
						/* Second invariant of viscoplastic strain calc, check */
						e1=Marke[mm1]-sbrit/2.0/markgg[mm2];
						/* Correction of divergence rate for plastic strain rate */
						if(e1<markf1[mm2])
						{
							e1=epsin-sbrit/2.0/nuduct;
							if(e1) *mdi0=2.0*e1*markf0[mm2]*timestepe;
						}
					}
					
					/* ! Recompute stress ! So stress no longer exceed strength */
					if(sigin && sbrit<sigin)
					{
						Markxx[mm1] *= sbrit/sigin;
						Markxy[mm1] *= sbrit/sigin;
						sigin=sbrit;
					}
					
					/* ! Recompute viscosity ! So decrease viscosity accordingly to localize deformation */
					nubrit=sbrit/(2.0*epsin+(sigin-sbrit)/timestepe/markgg[mm2]);
					if(nubrit<nueff) nueff=nubrit;
					
					/* Set initial plastic strain */
					if(Marke[mm1]<=0) Marke[mm1]=1e-20;
				}
			}
			else
			{
				if(yn==1) Marke[mm1]=0;
			}
		}
	}

	/* ------------------ End calculation viscoplastic viscosity ---------------------------- */
	
	/* Check calculated viscosity to be within hard code minimum and maximum */
	if(nueff<nubeg) nueff=nubeg; if(nueff>nuend) nueff=nuend;
	if(nueff<markn0[mm2]) nueff=markn0[mm2]; if(nueff>markn1[mm2]) nueff=markn1[mm2];
	
	// Pass final viscosity back to main model
	*mnu = nueff;
}
/* End OMP: Nu calc after reological equation */


/* Number of nearest left vertical line find */
long int m1serch(double cmx)
	/* cmx - X coordinate */
{
	/* Variables */
	long int m1,m10=0,m11=xnumx-1;

	/* Serch cycle */
	do
	{
		m1=(m10+m11)/2;
		if (gx[m1]>cmx) m11=m1; else m10=m1;
	}
	while((m11-m10)>1);
	if(m10>xnumx-2) m10=xnumx-2;

	return m10;
}
/* Number of nearest left vertical line find */


/* Number of nearest upper horizontal line find */
long int m2serch(double cmy)
	/* cmy - Y coordinate */
{
	/* Variables */
	long int m2,m20=0,m21=ynumy-1;

	/* Serch cycle */
	do
	{
		m2=(m20+m21)/2;
		if (gy[m2]>cmy) m21=m2; else m20=m2;
	}
	while((m21-m20)>1);
	if(m20>ynumy-2) m20=ynumy-2;

	return m20;
}
/* Number of nearest upper horizontal line find */


/* Erosion/Sedimentation Function for markers */
/* mardy - marker vertical size, m */
void erosmarkomp(long int mm1, int yn, long int m10, double x, double y, char markt[], double marke[], double markd[])
	/* mm1 - marker number */
	/* yn - current sedimnts type 2,3 */
	/* m1 - Up Left Node X,Y Num */
{
	/* Variables */
	double e,e0;

	/* Surface level elevation definition */
	/* Relativ Normalized coord Calc */
	e=(x-gx[m10])/(gx[m10+1]-gx[m10]);

	/* Surface level elevation for marker definition */
	e0=(e*ep[m10+1]+(1.0-e)*ep[m10]);

	/* Marker surface elevation definition */
	if(markt[mm1]<2)
	{
		/* Water/Air -> Sediments conversion */
		if(y>e0) {markt[mm1]=yn; marke[mm1]=0; markd[mm1]=-1.0;}        
	}
	if(markt[mm1]>1)
	{
		/* Rock->Water/Air conversion */
		if(y<e0) {markt[mm1]=0; marke[mm1]=0; markd[mm1]=-1.0;}
	}
}
/* OMP End Erosion/Sedimentation Function for markers */


/* OMP Rock to rock+melt transformation */
void meltingomp(double mtk, double mpb, long int mm1, int mm2, char Markt[], double Marke[], double *mxmelt, double *mhlatent)
	/* mtk - T, K */
	/* mpb - P, bar */
	/* mm1 - mark number */
{
	
	/* Melting related cahnge of the marker type */
	/* Check marker type */
	if (mm2==3 || mm2==4 || mm2==5 || mm2==6 || mm2==7 || mm2==8 || mm2==11 || mm2==16 || mm2==23 || mm2==24 || mm2==25 || mm2==26 || mm2==27 || mm2==28 || mm2==34 || mm2==36 || mm2==37 || mm2==38)
		if (mpb<0) mpb=0;
	switch(mm2)
	{
		/* Sediments, upper crust */
		case  3:
		case  4:
		case  5:
		case  17:
		case  23:
		case  24:
		case  25:
		case  26:
		case  37:
		/* Basalt, Gabbro */
		case  7:
		case  8:
		case  16:
		case  6:
		case  18:
		case  27:
		case  28:
		case  36:
		case  38:
			// mxmelt and mhlatent are already pointers to mem address, so you can enter them without &
			meltpart1omp(mtk,mpb,mm2,mxmelt,mhlatent);
			if(*mxmelt>0 && mm2<20) {Markt[mm1]+=20; Marke[mm1]=0;}
			if(*mxmelt<=0 && mm2>20) {Markt[mm1]-=20; Marke[mm1]=0;}
			return;

		/* Hydrated Peridotite */
		case  11:
		case  34:
			meltpart1omp(mtk,mpb,mm2,mxmelt,mhlatent);
			if(*mxmelt>0 && mm2==11) {Markt[mm1]=34; Marke[mm1]=0;}
			if(*mxmelt<=0 && mm2==34) {Markt[mm1]=14; Marke[mm1]=0;}	
			return;
		/* Others */
		default: return;
	}
}
/* OMP End Rock to rock+melt transformation */


/* Melt fraction, density, viscosity, heat capacity calculation */
void meltpartomp(double mtk, double mpb, double x, double y, long int mm1, int mm2, double *mro,double *mbb, double *maa, double *mnu, double *mcp, double *mkt, double *mgg, double *mxmelt, double *mhlatent)
	/* mtk - T, K */
	/* mpb - P, bar */
	/* x,y - XY location of point for Vx,Vy calc */
	/* mm1 - mark number */
	/* mm2 - mark type */
{
	/* Val buffer */
	double xmelt=0,ival,dmpb,dmtk,sduct,nueff,smin,smax,nmin,nmax,cpadd=0,vx0,vy0,pr0,sp0,ee0;
	long int m1,m10,m20;
	double Mnu,mdi0;
	double p_pl_out,p_ga_in,rokf,p_sp_in,p_ol_out,p_pv_in,p_sp_out,p_st_in;
	m10=m1serch(x);
	
	/* Check marker type */
	if (mm2==23 || mm2==24 || mm2==25 || mm2==26 || mm2==27 || mm2==28 || mm2==34 || mm2==36 || mm2==37 || mm2==38)
	{
		/* Calculate melt fraction */
		// mxmelt and mhlatent are already pointers to mem address, so you can enter them without &
		meltpart1omp(mtk,mpb,mm2,mxmelt,mhlatent);
		xmelt = *mxmelt;

		/* Standard adiabatic term: al=bro/(1+bro*(Tk-298.15)) */
		*mbb=(markbb[mm2]*xmelt+markbb[mm2-20]*(1.0-xmelt))/(1.0-(markbb[mm2]*xmelt+markbb[mm2-20]*(1.0-xmelt))*(mtk-298.15));
		*maa=(markaa[mm2]*xmelt+markaa[mm2-20]*(1.0-xmelt))/(1.0+(markaa[mm2]*xmelt+markaa[mm2-20]*(1.0-xmelt))*(mpb-1.0)*1e-3);

		/* Density */
		/* Ro=ro0 */
		if (densimod==0) 
		{
			*mro=markro[mm2]*xmelt+markro[mm2-20]*(1.0-xmelt);
		}
		/* Ro=ro0*(1-bro*(TK-298.15))*(1+aro*(Pkbar-0.001)) */
		else
		{
			ival=1.0;
	
	            /* ========================= */
                    /* Mantle phase transitions  */
                    /* ========================= */
			/*
			if(mm2>=29 && mm2<=34 && markex[mm1]>0) ival=1.0-0.04*markex[mm1];
			*/
			/* Eclogitization, St, Pv transitions in oceanic crust */
			if(mm2>=27 && mm2<=28)
				{
				/* Eclogitization Ito and Kennedy, 1971 */
	        		/*basalt=>garnet granulite (Ga-In) transition*/
	        		p_ga_in=-9222.0+mtk*14.0;
				/*Not to have granulites at pressure lower than 2 kbar*/
	       			if(p_ga_in<2000.0) p_ga_in=2000.0;
	     			/*garnet granulite=>eclogite (Pl-Out) transition*/
	        		p_pl_out=-1460.0+mtk*20.0;
	      			/*Not to have eclogites at pressure lower than 12 kbar*/
	     			if(p_pl_out<12000.0) p_pl_out=12000.0;
	 			if(mpb>p_ga_in)
				       	{
					rokf=0;
					if(mtk>teclmin)
						{
						if(mtk>teclmax)
	                            			{
	                            			rokf=0.16;
	                            			}
						else
	                            			{
	                            			rokf=0.16*(mtk-teclmin)/(teclmax-teclmin);
	                            			}
	                			}
				        if(mpb>=p_pl_out)
	                   			{
				               	ival=1.0+rokf;
	                       			}
	                   		else
	                        		{
	                        		ival=(1.0+rokf*(mpb-p_ga_in)/(p_pl_out-p_ga_in));
	                        		}
	                		}
				/* Coe->St transition Gerya et al., 2004, PCM */
	        		p_st_in=59100.0+mtk*22.6;
	        		if(mpb>p_st_in) ival*=1.06;
				/* Pv transition, Mishin et al., 2008 with slope from Ito et al., 1990 */
	        		/* Sp-out transition*/
	        		p_sp_out=354000.0-mtk*40.0;
	        		/* Pv-in transition*/
	        		p_pv_in=352000.0-mtk*40.0;
	        		if(mpb>p_pv_in)
	                		{
					rokf=0.08;
	                		if(mpb>=p_sp_out)
	                   			{
	                       			ival*=1.0+rokf;
	                        		}
	                		else
	                        		{
	                        		ival*=(1.0+rokf*(mpb-p_pv_in)/(p_sp_out-p_pv_in));
	                        		}
	                		}
	        		}
			/* Ol-Sp and Pv transitions in the mantle */
			if(mm2>=29 && mm2<=34) 
	        		{
				/* Ol-Sp transition, Katsura & Ito, 1989 */
	        		/* Ol-out transition*/
	        		p_ol_out=91000.0+mtk*27.0;
	        		/* Sp-in transition*/
	        		p_sp_in=66000.0+mtk*39.0;
	        		/*Limit width of Sp-Ol transition to 2 kbar */
	        		if(p_sp_in>p_ol_out-2000.0) p_sp_in=p_ol_out-2000.0;
	        		if(mpb>p_sp_in)
	                		{
					rokf=0.06;
	                		if(mpb>=p_ol_out)
	                   			{
	                       		ival*=1.0+rokf;
	                        		}
	                		else
	                        		{
	                        		ival*=(1.0+rokf*(mpb-p_sp_in)/(p_ol_out-p_sp_in));
	                        		}
	                		}
				/* Pv transition, Ito et al., 1990 */
	        		/* Sp-out transition*/
	        		p_sp_out=304000.0-mtk*40.0;
	        		/* Pv-in transition*/
	        		p_pv_in=302000.0-mtk*40.0;
	        		if(mpb>p_pv_in)
	                		{
					rokf=0.11;
	                		if(mpb>=p_sp_out)
	                   			{
	                       			ival*=1.0+rokf;
	                        		}
	                		else
	                        		{
	                        		ival*=(1.0+rokf*(mpb-p_pv_in)/(p_sp_out-p_pv_in));
	                        		}
	                		}
				}
			/* Density calculation with corrections */
			*mro=xmelt*markro[mm2]*(1.0-markbb[mm2]*(mtk-298.15))*(1.0+markaa[mm2]*(mpb-1.0)*1e-3)+(1.0-xmelt)*ival*markro[mm2-20]*(1.0-markbb[mm2-20]*(mtk-298.15))*(1.0+markaa[mm2-20]*(mpb-1.0)*1e-3);
			}
		/**/	
		/* Viscosity */
		/* Effective NU calc check */
		/* Little melt */
		// Assume similar to no melt, since go into viscalc..
		if(xmelt<0.1)
		{
			// QUESTION TARAS - why plastic reset here? (i switched yn=1 to yes wrt old version, but before was set to 0 here)
			// while mm2 going in is mm2-20 ? And mm2>20 returns immediately; ok that put here mm2-20 ?
	                viscalcomp(mtk,mpb,markx[mm1],marky[mm1],markv[mm1],markwa[mm1],markk[mm1],markp[mm1],markt[mm1],markexx[mm1],markexy[mm1],markxx,markxy,marke,mm1,mm2-20,1,m10,&Mnu,&mdi0);
                        *mnu=Mnu;
                        *mgg=markgg[mm2-20];
		}

		/* Significant melt */
		// Allowed to drop viscosity below minimum for rock type (init.t3c), but not below minimum for whole model (mode.t3c)
		else
		{
			/* Set viscosity and stress limits */
			nmin=MAXV(markn0[mm2],nubeg);
			nmax=MINV(markn1[mm2],nuend);
			smin=MAXV(marks0[mm2],strmin);
			smax=MINV(marks1[mm2],strmax);
			/* Calc effective strain rate after second strain rate tensor invariant EEii=(1/2SUM(EPSik^2))^(1/2) */
			m20=m2serch(y);
			allinteriomp(x,y,m10,m20,&vx0,&vy0,&pr0,&sp0,&ee0);
			// ee0=pow(eps[6]*eps[6]+eps[4]*eps[4],0.5); (was epsin)
			/* Effective NU calc check */
			nueff=marknu[mm2]*exp(2.5+pow((1.0-xmelt)/xmelt,0.48)*(1.0-xmelt));
			if(nueff<nmin) nueff=nmin;
			if(nueff>nmax) nueff=nmax;
			/* Ductile stress calc check */
			sduct=nueff*2.0*ee0;
			if(sduct<smin && ee0) {nueff=0.5*smin/ee0; sduct=smin;}
			if(sduct>smax) {nueff=0.5*smax/ee0; sduct=smax;}
			*mnu=nueff;
			/* Shear modulus */
			*mgg=markgg[mm2];
		}

		/* Heat capacity */
		*mcp=markcp[mm2]*xmelt+markcp[mm2-20]*(1.0-xmelt);

		/* heat conductivity */
		*mkt=((markkt[mm2]+markkf[mm2]/(mtk+77.0))*exp(markkp[mm2]*mpb))*xmelt+((markkt[mm2-20]+markkf[mm2-20]/(mtk+77.0))*exp(markkp[mm2-20]*mpb))*(1.0-xmelt);

		/* Additional melting adiabatic term, heat capacity */
		if(xmelt>0 && xmelt<1.0)
		{
			/* Melting adiabatic term: alm=-ro*(dHlat/dP)/T */
			/* Numerical differentiation */
			dmpb=mpb*0.001;
			meltpart1omp(mtk,mpb-dmpb,mm2,mxmelt,mhlatent);
			ival= *mhlatent;
			meltpart1omp(mtk,mpb+dmpb,mm2,mxmelt,mhlatent);
			ival-= *mhlatent;
			ival *= *mro / (mtk*2.0*dmpb*1e+5);
			*mbb+=ival;

			/* Melting heat capacity term: cpm=dHlat/dT */
			/* Numerical differentiation */
			dmtk=1.0;
			meltpart1omp(mtk+dmtk,mpb,mm2,mxmelt,mhlatent);
			ival= *mhlatent;
			meltpart1omp(mtk-dmtk,mpb,mm2,mxmelt,mhlatent);
			ival-= *mhlatent;
			ival/=2.0*dmtk;
			*mcp+=ival;
		}
	}
	else
	{
		*maa= *mbb= *mxmelt= *mhlatent= *mro= *mnu= *mcp= *mkt= 0;
	}
}
/* End OMP Rock to rock+melt transformation */




/* Melt fraction, latent heat calculation */
void meltpart1omp(double mtk, double mpb, int mm2, double *mxmelt, double *mhlatent)
	/* mtk - T, K */
	/* mpb - P, bar */
	/* x,y - XY location of point for Vx,Vy calc */
	/* mm1 - mark number */
	/* mm2 - mark type */
	/* yn  - type of calculation: 0 - Ro, 1 - Nu, 2 - Cp, 3 - kt */
{
	/* Val buffer */
	double xmelt=0,hlatent=0,ival;
	long int m1;
	double ykm=mpb*3e-3,ts=0,tl=0;
	
	/* Calculate melt fraction using marker type */
	if (ykm>0)
		switch(mm2)
	{
			/* Sediments: latent heat 300 kJ/kg (Bittner & Schmeling, 1995) */
		case  3:
		case  4:
		case  5:
		case 17:
		case 23:
		case 24:
		case 25:
		case 37:
			/* Wet Solidus Temperature, Johannes, 1985, Poli & Schmidt, 2002 */
			if (ykm<36.0) 
			{
				ts=889.0+536.6/(ykm+1.609)+18.21/(ykm+1.609)/(ykm+1.609);
			}
			else
			{
				ts=831.3+2.0*ykm;
			}
			/* Dry Granite Liquidus, Johannes, 1985 */
			tl=1262.0+3.0*ykm;
			hlatent=300000.0;
			break;

			/* Basalt, Gabbro: latent heat 380 kJ/kg (Bittner & Schmeling, 1995) */
		case 7:
		case 8:
		case 16:
		case 27:
		case 28: 
		case 36: 
		case  6:
		case 18:
		case 26:
		case 38:
			/* Wet solidus, Schmidt & Poli, 1998  */
			if (ykm<48.0) 
			{
				ts=972.6-2111.0/(ykm+10.63)+70033.0/(ykm+10.63)/(ykm+10.63);
			}
			else
			{
				ts=935.4+0.1162*ykm+0.006937*ykm*ykm;
			}
			/* Dry Toleitic Basalt Liquidus, Hess, 1989 */
			tl=1423.15+3.5*ykm;
			hlatent=380000.0;
			break;

		/* Peridotite: latent heat 400 kJ/kg Turcotte & Schubert, 1982, p.171 */
		case 11:
		case 34:
			/* Wet solidus, Schmidt & Poli, 1998  */
			if (ykm<72.0) 
			{
				ts=1239.8+1493.0/(ykm+9.701);
			}
			else
			{
				ts=1266.3-0.3948*ykm+0.003893*ykm*ykm;
			}
			/* Dry Peridotite Liquidus, Hess, 1989 */
			tl=2073.15+3.8*ykm;
			hlatent=400000.0;
			break;

		/* Other rocks - No melting */
		default:
		break;
	}

	/* Melt fraction, latent heat calculation */
	*mxmelt = *mhlatent = 0;
	if(tl)
	{
		/* Melt fraction calc, check */
		xmelt=(mtk-ts)/(tl-ts);
		if(xmelt<0) xmelt=0;
		if(xmelt>1.0) xmelt=1.0;
		*mxmelt = xmelt;
		
		/* Latent heat calc */
		hlatent *= xmelt;
		*mhlatent=hlatent;
	}
}
/* End OMP Melt fraction, latent heat calculation */


/* Hydration front progress after H2O budget */
double hydration2omp()
{
	/* Val buffer */
	double ysurf,vfiltr,yfiltr,dydx,dydx1,sy1,sy2,sy3,sy4,sy5,e1,mwamin,x0,y0,x1,y1,vx1,vy1;
	double hytimesum,hytimesum0;
	/* TD Database variables */
	double W0,W1,W2,W3,R0,R1,R2,R3,n,e,dx,dy;
	double mtk,mpb,mwa,mro,dmwa,wro;
	double Mgg,Mro,Mwa,Mcp,Mbb,Maa,Mdhh,Mkt;
	long int m1,m2,m3,mm1,marknum1=marknum;
	int mm2,mm3,n1,n2;

	fprintf(fp_log,"\n WATER Transport BEGIN \n");fflush(fp_log);
	
	/* Marker steps */
	dx=dxwater;
	dy=dywater;

	/* Min water contents in the hydraten mantle wt% */
	mwamin=0.1;
	/* Min Distance from erosion surface for water release */
	ysurf=8000.0;
	
	/* Clear wa[] wt */
	for (m1=0;m1<nodenum;m1++)
	{
		wa0[m1]=0;
		wa1[m1]=0;
		sol0[m1]=0;
		sol1[m1]=0;
		sol0[nodenum+m1]=1e+30;
		sol1[nodenum+m1]=-1e+30;
		sol0[nodenum2+m1]=1e+30;
		sol1[nodenum2+m1]=-1e+30;
		fre0[         m1]=1e+30;
		fre0[nodenum +m1]=-1e+30;
		fre0[nodenum2+m1]=1e+30;
		fre0[nodenum3+m1]=-1e+30;
	}
	
	/* Fluid marker generation cycle */
	double start=omp_get_wtime();
	for (mm1=0;mm1<marknum;mm1++)
	{
		// Reset fluid presence indicator for next marker for loop
		markwa[mm1] = 0;

		/* Marker type */
		mm2=(int)markt[mm1]; if (mm2>=100) mm2-=100;
		
		/* Marker cell number */
		m1=m1serch(markx[mm1]);
		m2=m2serch(marky[mm1]);
		m3=m1*ynumy+m2;
		
		/* Erosion surface */
		e1=(markx[mm1]-gx[m1])/(gx[m1+1]-gx[m1]);
		sy1=(e1*ep[m1+1]+(1.0-e1)*ep[m1]);

		/* Check markers out of grid and within hydration range */
		if(markx[mm1]>0 && marky[mm1]>0 && (markx[mm1])<xsize && (marky[mm1])<ysize && (markk[mm1]>0 || markt[mm1]>=50) && markt[mm1]<100)
			if((markd[mm1])>=0 && (markw[mm1])>=0 && mm2>1 && mm2!=9 && mm2!=10)
		{
			if(mm2<50)
			{
				/* P, T parameters calc */
				mpb=1e-5*allinterpomp(markx[mm1],marky[mm1],m1,m2);
				mtk=(markk[mm1]);

				/* Mantle to Antigorite transformation */
				antigoromp(mtk,mpb,markx[mm1],marky[mm1],mm1,m1,markt);

				/* Rocks to rock+melt transformation */
				if (markt[mm1]>=20)
				{
					/* Check melting extent */
					if(fre0[        +m3]>markx[mm1]-dx) fre0[         m3]=markx[mm1]-dx;
					if(fre0[nodenum +m3]<markx[mm1]+dx) fre0[nodenum +m3]=markx[mm1]+dx;
					if(fre0[nodenum2+m3]>marky[mm1]-dy) fre0[nodenum2+m3]=marky[mm1]-dy;
					if(fre0[nodenum3+m3]<marky[mm1]+dy) fre0[nodenum3+m3]=marky[mm1]+dy;
				}
				
				/* Compute TD variables */
				tdbasecalcomp(markx[mm1],marky[mm1],mtk,mpb,mm2,mm1,m1,&Mgg,&Mro,&Mwa,&Mcp,&Mbb,&Maa,&Mdhh,&Mkt);
				mro=Mro;
				mwa=Mwa;

				/* Water changes in kg/m3 calc */
				dmwa=mro*(mwa-markw[mm1])*1e-2;
				//{fprintf(fp_log,"H2O MARKER %ld %d %d %e %e   %e %e  %e %e   %e",mm1,mm2,mm3,mtk-273.15,mpb/1000.0,mwa,mro,markw[mm1],markd[mm1],dmwa);getchar();}
				//{fprintf(fp_log,"H2O RELEASE %ld %d %d %e %e   %e %e  %e %e   %e",mm1,mm2,mm3,mtk-273.15,mpb/1000.0,mwa,mro,markw[mm1],markd[mm1],dmwa);getchar();}

				/* Add water changes to the current cell, kg/m3 */
				/* Water release */
				if ((markw[mm1]-mwa)>dmwamin)
				{
					/* Save new water content */
					markw[mm1]=mwa;
					/* Generation of fluid marker (NO FLUID From melts */
					if (markt[mm1]<20 && marky[mm1]>sy1)
					{
						markt[marknum1]=markt[mm1]+50;
						markx[marknum1]=markx[mm1];
						marky[marknum1]=marky[mm1];
						markk[marknum1]=markk[mm1];
						markd[marknum1]=1050.0;
						markw[marknum1]=-dmwa;
						/* Add aditional markers counter */
						marknum1++;
							
						// If new marker is interesting for picking algorithm, flag to follow
						// Note is hard-coded in i2.c as well. Only here excluded fluid markers, since immobile can not become fluid
						if ( start_cond==1 && marky[marknum1]<85e3 && markx[marknum1]>gx[m10_hr] && markx[marknum1]<gx[m11_hr] && markt[marknum1]>49 && markt[marknum1]<100)
						{
							follow[marknum1]=2;
							nmf++;
						}
							
						/* Check hydration extent */
						if(sol0[nodenum+m3]>markx[mm1]-dx) sol0[nodenum+m3]=markx[mm1]-dx;
						if(sol1[nodenum+m3]<markx[mm1]+dx) sol1[nodenum+m3]=markx[mm1]+dx;
						if(sol0[nodenum2+m3]>marky[mm1]-dy) sol0[nodenum2+m3]=marky[mm1]-dy;
						if(sol1[nodenum2+m3]<marky[mm1]+dy) sol1[nodenum2+m3]=marky[mm1]+dy;
					}
				}
				else
					/* Water consuming */
				{
					if(dmwa>0)
					{
						wa1[m3]+=dmwa;
						sol1[m3]+=1.0;
					}
				}
			}
			else
				/* Fluid marker count */
			{
				/* Check position */
				if(marky[mm1]>sy1)
				{
					/* Check hydration extent */
					if(sol0[nodenum+m3]>markx[mm1]-dx) sol0[nodenum+m3]=markx[mm1]-dx;
					if(sol1[nodenum+m3]<markx[mm1]+dx) sol1[nodenum+m3]=markx[mm1]+dx;
					if(sol0[nodenum2+m3]>marky[mm1]-dy) sol0[nodenum2+m3]=marky[mm1]-dy;
					if(sol1[nodenum2+m3]<marky[mm1]+dy) sol1[nodenum2+m3]=marky[mm1]+dy;
				}
				else
					/* Erase fluid marker */
				{
					markx[mm1]=-1.0;
					markk[mm1]=0;
				}
			}
		}
	}

	/* Rock hydration cycle: rocks get hydrated by changing marker type mm2 */
	start=omp_get_wtime();
	for (mm1=0;mm1<marknum;mm1++)
		if(markx[mm1]>0 && marky[mm1]>0 && (markx[mm1])<xsize && (marky[mm1])<ysize && markt[mm1]<50)
	{
		/* Marker cell number */
		m1=m1serch(markx[mm1]);
		m2=m2serch(marky[mm1]);
		m3=m1*ynumy+m2;
		
		/* Check markers within hydration range */
		if(markx[mm1]>sol0[nodenum+m3] && marky[mm1]>sol0[nodenum2+m3] && (markx[mm1])<sol1[nodenum+m3] && (marky[mm1])<sol1[nodenum2+m3])
		{
			/* Fluid presence mark */
			markwa[mm1]=1;
			if(markt[mm1]==9 || markt[mm1]==10 || markt[mm1]==12 || markt[mm1]==14 || markt[mm1]==5 || markt[mm1]==6)
			{
				/* Mantle Hydration */
				if (markt[mm1]!=5 && markt[mm1]!=6)
				{
					mm2=markt[mm1]=11;
				}
				else
				{
					mm2=markt[mm1]=markt[mm1]+12;
				}
				/* P, T parameters calc */
				mpb=1e-5*allinterpomp(markx[mm1],marky[mm1],m1,m2);
				mtk=(markk[mm1]);

				/* Mantle to Antigorite transformation */
				antigoromp(mtk,mpb,markx[mm1],marky[mm1],mm1,m1,markt);

				/* Rocks to rock+melt transformation */
				if (markt[mm1]>=20)
				{
					/* Check melting extent */
					if(fre0[        +m3]>markx[mm1]-dx) fre0[         m3]=markx[mm1]-dx;
					if(fre0[nodenum +m3]<markx[mm1]+dx) fre0[nodenum +m3]=markx[mm1]+dx;
					if(fre0[nodenum2+m3]>marky[mm1]-dy) fre0[nodenum2+m3]=marky[mm1]-dy;
					if(fre0[nodenum3+m3]<marky[mm1]+dy) fre0[nodenum3+m3]=marky[mm1]+dy;
				}

				/* Thermodynamic database use for Ro as function of Water content */
				/* Compute TD variables */
				tdbasecalcomp(markx[mm1],marky[mm1],mtk,mpb,mm2,mm1,m1,&Mgg,&Mro,&Mwa,&Mcp,&Mbb,&Maa,&Mdhh,&Mkt);
				mro=Mro;
				mwa=Mwa;

				/* Water changes in kg/m3 calc */
				dmwa=mro*(mwa-markw[mm1])*1e-2;

				/* Add water changes to the current cell, kg/m3 */
				/* Water consuming */
				if (dmwa>0)
				{
					wa1[m3]+=dmwa;
					sol1[m3]+=1.0;
				}
			}
		}
	}

	/* Fluid marker computing cycle */
	start=omp_get_wtime();
	for (mm1=0;mm1<marknum1;mm1++)
	{
		/* Check markers out of grid and within hydration range */
		if(markt[mm1]>=50 && markt[mm1]<100 && markx[mm1]>0 && marky[mm1]>0 && (markx[mm1])<xsize && (marky[mm1])<ysize)
		{
			/* Marker cell number */
			m1=m1serch(markx[mm1]);
			m2=m2serch(marky[mm1]);
			m3=m1*ynumy+m2;

			/* Erosion surface */
			e1=(markx[mm1]-gx[m1])/(gx[m1+1]-gx[m1]);
			sy1=(e1*ep[m1+1]+(1.0-e1)*ep[m1]);
			/* Water in melt region conversion */
			if(markd[mm1]<1100.0 && markx[mm1]>fre0[m3] && marky[mm1]>fre0[nodenum2+m3] && markx[mm1]<fre0[nodenum+m3] && marky[mm1]<fre0[nodenum3+m3]) markd[mm1]=1150.0;

			/* Check position, no fluid above erosion/sedimentation level, no fluid passing through the melt */
			if(marky[mm1]>sy1 && marky[mm1]<zdeep && (markd[mm1]<1100.0 || (markx[mm1]>fre0[m3] && marky[mm1]>fre0[nodenum2+m3] && markx[mm1]<fre0[nodenum+m3] && marky[mm1]<fre0[nodenum3+m3])))
			{
				wa0[m3]+=markw[mm1];
				sol0[m3]+=1.0;
			}
			else
				/* Erase fluid marker */
			{
				markx[mm1]=-1.0;
				markk[mm1]=0;
			}
		}
	}
	if (printmod==10000) fprintf(fp_log,"\n  Time taken for fluid computing cycle = %e s \n",omp_get_wtime()-start);

	/* Fluid marker consuming cycle */
	start=omp_get_wtime();
	for (mm1=0;mm1<marknum1;mm1++)
	{
		/* Marker type */
		mm2=(int)markt[mm1]; if (mm2>=100) mm2-=100;
		// What use? since will not use mm1>100 anyway..

		/* Marker cell number */
		m1=m1serch(markx[mm1]);
		m2=m2serch(marky[mm1]);
		m3=m1*ynumy+m2;

		/* Change water consuming rocks  and fluid makers */
		if(markx[mm1]>0 && marky[mm1]>0 && (markx[mm1])<xsize && (marky[mm1])<ysize && (markk[mm1]>0 || markt[mm1]>=50) && markt[mm1]<100)
			if((markd[mm1])>=0 && (markw[mm1])>=0 && mm2>1 && mm2!=9 && mm2!=10 && mm2!=12 && mm2!=14 && mm2!=5 && mm2!=6)
		{
			// For all assimilating rock types: 0-50, except those one line above
			if(mm2<50)
			{
				/* P, T parameters calc */
				// Why need to do this every time again?
				mpb=1e-5*allinterpomp(markx[mm1],marky[mm1],m1,m2);
				mtk=markk[mm1];

				/* Thermodynamic database use for Ro, Water */
				/* Compute TD variables */			
				tdbasecalcomp(markx[mm1],marky[mm1],mtk,mpb,mm2,mm1,m1,&Mgg,&Mro,&Mwa,&Mcp,&Mbb,&Maa,&Mdhh,&Mkt);
				mwa=Mwa;
				mro=Mro;

				/* Water change */
				dmwa=mwa-markw[mm1];

				/* Add water changes to the current cell, kg/m3 */
				/* Water consuming */
				if(dmwa>0)
				{
					if (wa1[m3]<=wa0[m3])
					{
						/* Save complete new water content */
						markw[mm1]=mwa;
					}
					else
					{
						/* COmpute, Save partial new water content */
						markw[mm1]=markw[mm1]+dmwa*wa0[m3]/wa1[m3];
					}
				}
			}
			// For all fluid markers: 50-100
			else
				/* Fluid marker change */
			{
				// Evaluate wether all free water is finished 
				if(wa1[m3]<wa0[m3])
				{
					/* Count water changes for fluid marker */
					markw[mm1]*=1.0-wa1[m3]/wa0[m3];
				}
				else
					/* Erase fluid marker */
				{
					markx[mm1]=-1.0;
					markk[mm1]=0;
				}
			}
		}
	}

	/* Reset aditional markers */
	fprintf(fp_log,"\n WATER BEG Number of markers: OLD = %ld NEW = %ld \n",marknum,marknum1); fflush(fp_log);
	mm1=0;
	while(marknum1>marknum && mm1<marknum)
	{
		/* Reload marker */
		if((markx[mm1]<0 || marky[mm1]<0 || (markx[mm1])>xsize || (marky[mm1])>ysize) && markt[mm1]<100) 
		{
			/* Decrease aditional markers counter */
			marknum1--;
			if(markx[marknum1]>=0);
			{
				/* Type save */
				markt[mm1]=markt[marknum1];
				/* X,Y, water reload */
				markx[mm1]=markx[marknum1];
				marky[mm1]=marky[marknum1];
				markw[mm1]=markw[marknum1];
				markd[mm1]=markd[marknum1];
				markk[mm1]=markk[marknum1];
			}
		}
		/* Increase markers counter */
		mm1++;
	}
	fprintf(fp_log,"\n WATER END Number of markers: OLD = %ld NEW = %ld \n",marknum,marknum1); fflush(fp_log);
	/* Set new marker number */
	marknum=marknum1;

	return 0;
}
/* End OMP Hydration front progress after H2O budget */


/* Erosion Surface progress */
void erosion()
{
	/* Val buffer */
	double v0,v1,dydx,x1,vx1,vy1,dy;
	double ertimesum,ertimesum0;
	long int m1,m2;
	/**/
	/* Erosion Solution Cycle ------------------------------------------ */
	ertimesum=0;
	ertimesum0=timestep;
	do
	{
		/* Save old cycle results */
		for (m1=0;m1<xnumx;m1++)
		{
			ep0[m1]=ep[m1];
			ep0[xnumx+m1]=ep[xnumx+m1];
		}
		/**/
		/**/
		/**/
		/* Initial timestep definition */
		timestep=ertimesum0-ertimesum;
		/**/
		/**/
		/**/
		/* Erosion timestep definition using material velosity field */
		for (m1=0;m1<xnumx;m1++)
		{
			/* Calc horisontal Coordinate */
			x1=gx[m1];
			/**/
			/* EROSION SURFACE */
			/* Calc material velocity on the Surface using velosity field */
			allinteri(x1,ep0[m1]);
			vx1=eps[11];
			vy1=eps[12];
			/* Check horizontal timestep */
			/* Calc x derivative of y position of the Surface using upwind differences */
			dydx=0;
			if(vx1>0 && m1>0)
			{
				timestep=MINV(timestep,(gx[m1]-gx[m1-1])/vx1);
				/*
				fprintf(fp_log,"111 %ld %e %e %e %e %e %e %e",m1,vx1,vy1,(gx[m1]-gx[m1-1])/vx1,ertimesum0,ertimesum,ertimesum0-ertimesum,timestep);getchar();
				*/
			}
			if(vx1<0 && m1<xnumx-1)
			{
				timestep=MINV(timestep,(gx[m1]-gx[m1+1])/vx1);
				/*
				fprintf(fp_log,"222 %ld %e %e %e %e %e %e %e",m1,vx1,vy1,(gx[m1]-gx[m1+1])/vx1,ertimesum0,ertimesum,ertimesum0-ertimesum,timestep);getchar();
				*/
			}
			/* Check vertical timestep */
			if(vy1)
			{
				/* Horizontal line num definition */
				m2=m2serch(ep0[m1]);
				/* Check timestep */
				timestep=MINV(timestep,(gy[m2+1]-gy[m2])/ABSV(vy1));
				/*
				fprintf(fp_log,"333 %ld  %e %e %e %e %e %e %e",m2,vx1,vy1,(gy[m2+1]-gy[m2])/ABSV(vy1),ertimesum0,ertimesum,ertimesum0-ertimesum,timestep);getchar();
				*/
			}
			/**/
			/**/
			/* INITIAL SURFACE */
			/* Calc material velocity on the Initial Surface using velosity field */
			allinteri(x1,ep0[xnumx+m1]);
			vx1=eps[11];
			vy1=eps[12];
			/* Check horizontal timestep */
			/* Calc x derivative of y position of the Surface using upwind differences */
			dydx=0;
			if(vx1>0 && m1>0)
			{
				timestep=MINV(timestep,(gx[m1]-gx[m1-1])/vx1);
				/*
				fprintf(fp_log,"444 %ld %e %e %e %e %e %e %e",m1,vx1,vy1,(gx[m1]-gx[m1-1])/vx1,ertimesum0,ertimesum,ertimesum0-ertimesum,timestep);getchar();
				*/
			}
			if(vx1<0 && m1<xnumx-1)
			{
				timestep=MINV(timestep,(gx[m1]-gx[m1+1])/vx1);
				/*
				fprintf(fp_log,"555 %ld %e %e %e %e %e %e %e",m1,vx1,vy1,(gx[m1]-gx[m1+1])/vx1,ertimesum0,ertimesum,ertimesum0-ertimesum,timestep);getchar();
				*/
			}
			/* Check vertical timestep */
			if(vy1)
			{
				/* Horizontal line num definition */
				m2=m2serch(ep0[xnumx+m1]);
				/* Check timestep */
				timestep=MINV(timestep,(gy[m2+1]-gy[m2])/ABSV(vy1));
				/*
				fprintf(fp_log,"666 %ld  %e %e %e %e %e %e %e",m2,vx1,vy1,(gy[m2+1]-gy[m2])/ABSV(vy1),ertimesum0,ertimesum,ertimesum0-ertimesum,timestep);getchar();
				*/
			}
		}
		/*
		fprintf(fp_log,"777 %e %e %e %e",ertimesum0,ertimesum,ertimesum0-ertimesum,timestep);getchar();
		*/
		/**/
		/**/
		/**/
		/* Displace Surface boundary */
		/*
		for (m1=1;m1<xnumx-1;m1++)
		*/
		for (m1=0;m1<xnumx;m1++)
		{
			/* EROSION SURFACE */
			/* Calculation of errosion rate */
			v0=0;
			if(ep0[m1]<eroslev)
			{
				v0=eroscon+eroskoe*(eroslev-ep0[m1]);
			}
			/* Calculation of sedimentation rate */
			v1=0;
			if(ep0[m1]>sedilev)
			{
				v1=sedicon+sedikoe*(ep0[m1]-sedilev);
			}
			/* Calc horisontal Coordinate */
			x1=gx[m1];
			/**/
			/* Calc material velocity on the Surface using velosity field */
			allinteri(x1,ep0[m1]);
			vx1=eps[11];
			vy1=eps[12];
			/**/
			/* Erase erosion/sedimentation rate for marginal points */
			if((m1==0 && vx1>0) || (m1==xnumx-1 && vx1<0)) v0=v1=0;
			/**/
			/* Calc x derivative of y position of the Surface using upwind differences */
			dydx=0;
			if(vx1>0 && m1>0)
			{
				dydx=(ep0[m1]-ep0[m1-1])/(gx[m1]-gx[m1-1]);
				/*
				fprintf(fp_log,"AAA %e %e",ep0[m1],dydx);getchar();
				*/
			}
			if(vx1<0 && m1<xnumx-1)
			{
				dydx=(ep0[m1+1]-ep0[m1])/(gx[m1+1]-gx[m1]);
				/*
				fprintf(fp_log,"BBB %e %e",ep0[m1],dydx);getchar();
				*/
			}
			/* Recalc new Surface position */
			ep[m1]+=timestep*(v0-v1+vy1-dydx*vx1);
			/*
			fprintf(fp_log,"SURFACE %ld %e %e %e %e %e %e %e %e",m1,x1,v0,v1,vx1,vy1,dydx,ep[m1]);getchar();
			*/
			/**/
			/**/
			/**/
			/* INITIAL SURFACE */
			/* Initial surface displacement */
			/* Calc material velocity on the Surface using velosity field */
			allinteri(x1,ep0[xnumx+m1]);
			vx1=eps[11];
			vy1=eps[12];
			/* Calc x derivative of y position of Initial Surface using upwind differences */
			dydx=0;
			if(vx1>0 && m1>0)
			{
				dydx=(ep0[xnumx+m1]-ep0[xnumx+m1-1])/(gx[m1]-gx[m1-1]);
				/*
				fprintf(fp_log,"AAA %e ",dydx);getchar();
				fprintf(fp_log,"AAA %e ",dydx);getchar();
				*/
			}
			if(vx1<0 && m1<xnumx-1)
			{
				dydx=(ep0[xnumx+m1+1]-ep0[xnumx+m1])/(gx[m1+1]-gx[m1]);
				/*
				fprintf(fp_log,"BBB %e ",dydx);getchar();
				*/
			}
			/* Recalc new Initial Surface position */
			ep[xnumx+m1]+=timestep*(vy1-dydx*vx1);
			/**/
		}
		/**/
		/**/
		/**/
		/**/
		/* Relax EROSION surface */
		if (0==0)
			for (m1=0;m1<xnumx-1;m1++)
		{
			/* Calc x derivative of y position */
			dydx=(ep[m1+1]-ep[m1])/(gx[m1+1]-gx[m1]);
			/* Relax surface for critical slope */
			if(dydx>slopemax)
			{
				dy=((ep[m1+1]-ep[m1])-slopemax*(gx[m1+1]-gx[m1]))/2.0;
				ep[m1]  +=dy;
				ep[m1+1]-=dy;
				/*
				dydx=(ep[m1+1]-ep[m1])/(gx[m1+1]-gx[m1]);
				fprintf(fp_log,"AAA %ld %e %e",m1,slopemax,dydx);getchar();
				*/
			}
			if(dydx<-slopemax)
			{
				dy=((ep[m1+1]-ep[m1])+slopemax*(gx[m1+1]-gx[m1]))/2.0;
				ep[m1]  +=dy;
				ep[m1+1]-=dy;
				/*
				dydx=(ep[m1+1]-ep[m1])/(gx[m1+1]-gx[m1]);
				fprintf(fp_log,"BBB %ld %e %e",m1,slopemax,dydx);getchar();
				*/
			}
		}
		/**/
		/**/
		/**/
		/* Add Erosion step */
		ertimesum+=timestep;
		/**/
		/**/
		/**/
		/* Print Results */
		if (printmod) { fprintf(fp_log,"\n EROSION STEP = %e YEARS    EROSION TIME = %e YEARS \n",timestep/3.15576e+7,ertimesum/3.15576e+7); fflush(fp_log); }
	}
	while(ertimesum<ertimesum0);
	/* Restore timestep */
	timestep=ertimesum0;
}
/* Erosion Surface progress */


/* Thermodynamic database use for ro, Cp */
// Within a loop over all markers, do: 
// Interpolation properties between four nearest points in thermodynamic database dep. on T,P,composition 
void tdbasecalcomp(double x, double y, double mtk, double mpb, int mm2, long int mm1, long int m10, double *Mgg, double *Mro, double *Mwa, double *Mcp, double *Mbb, double *Maa, double *Mdhh, double *Mkt)
{
	/* TD Database variables,  dTK,dPB - TK, PB step for tabulation in TD database */
	double H0,H1,H2,H3,R0,R1,R2,R3,G0,G1,G2,G3,W0,W1,W2,W3,n,e;
	/* Val Buffers */
	int n1,n2,mm3,ynpb;
	double mhh0,mhh1,mdhh,maa,mwa,dmwa,wro,mro,mcp,mbb,mgg,mkt,mkt1,pbmax,xold,kr01,kr1,kr10,xkr,krad;
	long int m1=m10;
	double sy1,e1;

	/* Maximal pressure for the shallow database */
	pbmax=pbmin+pbstp*(double)(pbnum-1);
	/* Adiabate computing */
	ynpb=0; if(1==0 && timesum<3.15576e+7*1e+3) {fprintf(fp_log,"in adiabate: can not right ? \n"); fflush(fp_log); mpb*=timesum/(3.15576e+7*1e+3); ynpb=1;}

	/* Reset TD variables */
	*Mgg=*Mro=*Mwa=*Mcp=*Mbb=*Maa=0;

	/* Thermal conductivity */
	/* m895 Dry peridotite Fe=12 */
	/* Olivine: Hoffmeister, 1999; Hoffmeister & Yuen, 2005 */
	if(mpb<235000.0)
	{
		/* Lattice k */
		mkt1=(1.878+770.9/MINV(mtk,1200.0))*(1.0+4.26e-6*mpb);
		/* Radiative k 0.1 mm */
		kr01=pow(mtk/4000.0,3.0);
		/* Radiative k 1 mm */
		kr1=pow(mtk/1774.0,3.0);
		/* Radiative k 10 mm */
		xkr=pow(mtk/1636.0,10.0);
		xkr/=xkr+1.0; kr10=pow((mtk-1000.0*xkr)/1011.0,3.0)-0.7713*xkr;
	}
	/* Perovskite: Hoffmeister, 1999; Hoffmeister & Yuen, 2005 */
	else
	{
		/* Lattice k */
		mkt1=(1.291+1157.0/MINV(mtk,2100.0))*(1.0+2.50e-6*mpb);
		/* Radiative k 0.1 mm */
		kr01=pow(mtk/3591.0,3.0);
		/* Radiative k 1 mm */
		kr1=pow(mtk/2117.0,3.0);
		/* Radiative k 10 mm */
		xkr=pow(mtk/1500.0,4.0); xkr/=xkr+1.0;
		kr10=pow((mtk+4000.0*xkr)/5776.0,3.0)+2.822*xkr;
	}
	krad=kr1;

	/* Shallow TD base type */
	if(mpb<pbmax && ynpb==0)
	{
		/* TD base type */
		switch (mm2)
		{
			/* Dry Upper crust */
			case 5: mm3=11; break;
			/* Wet Upper crust */
			case 17: mm3=12; break;
			/* Dry Lower crust */
			case 6: mm3=13; break;
			/* Wet Lower crust */
			case 18: mm3=14; break;
			/* Sediments */
			case 2:
			case 3:
			case 4: mm3=5; break;
			/* Molten Sediments */
			case 37:
			case 25:
			case 22:
			case 23:
			case 24: mm3=6; break;
			/* Basalt */
			case 16:
			case 7: mm3=7; break;
			/* Molten Basalt */
			case 36:
			case 27: mm3=8; break;
			/* Gabbro */
			case 38:
			case 26:
			case 8: mm3=3; break;
			/* Molten Gabbro */
			case 28: mm3=4; break;
			/* Dry peridotite */
			case 9:
			case 12:
			case 14:
			case 10: mm3=0; break;
			/* Wet peridotite */
			case 13:
			case 11: mm3=1; break;
			/* Molten peridotite */
			case 34: mm3=2; break;
			/* Unknown type */
			default: {fprintf(fp_log,"Shallow TD: Unknown rock type for TD database %d, for marker %ld with T= %f, P=%f \n",mm2,mm1,mtk,mpb); fflush(fp_log); exit(0);}
		}
		
		/* ABCD-4Cell Number */
		// Get weights for nearest points in thermodynamic database
		e=(mtk-tkmin)/tkstp;
		if(e<0) e=0;
		if(e>(double)(tknum-1)) e=(double)(tknum-1);
		n=(mpb-pbmin)/pbstp;
		if(n<0) n=0;
		if(n>(double)(pbnum-1)) n=(double)(pbnum-1);
		n1=(int)(e);
		if(n1>tknum-2) n1=tknum-2;
		n2=(int)(n);
		if(n2>pbnum-2) n2=pbnum-2;
		/* e,n Calc */
		e=(e-(double)(n1));
		n=(n-(double)(n2));
		/* Ro H values */
		/* 0 2 */
		/* 1 3 */
		R0=td[n1  ][n2  ][mm3][0]*1000.0;
		R1=td[n1  ][n2+1][mm3][0]*1000.0;
		R2=td[n1+1][n2  ][mm3][0]*1000.0;
		R3=td[n1+1][n2+1][mm3][0]*1000.0;
		H0=td[n1  ][n2  ][mm3][1]*1000.0*4.1837;
		H1=td[n1  ][n2+1][mm3][1]*1000.0*4.1837;
		H2=td[n1+1][n2  ][mm3][1]*1000.0*4.1837;
		H3=td[n1+1][n2+1][mm3][1]*1000.0*4.1837;
		W0=td[n1  ][n2  ][mm3][4];
		W1=td[n1  ][n2+1][mm3][4];
		W2=td[n1+1][n2  ][mm3][4];
		W3=td[n1+1][n2+1][mm3][4];
		G0=td[n1  ][n2  ][mm3][3]*1000.0;G0*=G0*R0;
		G1=td[n1  ][n2+1][mm3][3]*1000.0;G1*=G1*R1;
		G2=td[n1+1][n2  ][mm3][3]*1000.0;G2*=G2*R2;
		G3=td[n1+1][n2+1][mm3][3]*1000.0;G3*=G3*R3;
		/* Shear modulus calc by interpolation */
		mgg=((G0*(1.0-n)+G1*n)*(1.0-e)+(G2*(1.0-n)+G3*n)*e);
		/* Ro calc by interpolation */
		mro=((R0*(1.0-n)+R1*n)*(1.0-e)+(R2*(1.0-n)+R3*n)*e);
		/* Water wt% calc by interpolation */
		mwa=((W0*(1.0-n)+W1*n)*(1.0-e)+(W2*(1.0-n)+W3*n)*e);
		/* Add pore fluid */
		/* Erosion surface */
		e1=(x-gx[m10])/(gx[m10+1]-gx[m10]);
		sy1=y-(e1*ep[m10+1]+(1.0-e1)*ep[m10]);
		if(marks0[mm2]>0 && sy1>0 && sy1<zmpor && mtk<tkpor) 
		{
			dmwa=marks0[mm2]*(tkpor-mtk)/(tkpor-273.15)*(zmpor-sy1)/zmpor;
			mwa+=dmwa;
			wro=1050.0;
			mro=mro/(1.0+dmwa*1e-2*(mro/wro-1.0));
		}
		/* Cp calc by interpolation */
		mcp=((H2-H0)*(1.0-n)+(H3-H1)*n)/tkstp;
		if(mcp<1e+2) mcp=1e+2; else if(mcp>5e+4) mcp=5e+4;
		/* Effective adiabatic betta=1/V*dV/dT=ro/T*[-dH/dP+V] calc by interpolation */
		mbb=(2.0/(R1+R0)-(H1-H0)/pbstp/1e+5)*(1.0-e)+(2.0/(R3+R2)-(H3-H2)/pbstp/1e+5)*e;
		mbb*=mro/mtk;
		if(mbb<-1e-2) mbb=-1e-2; else if(mbb>1e-2) mbb=1e-2;
		/* Effective compressibility term alpha=1/ro*d(ro)/dP calc by interpolation */
		maa=(2.0/(R1+R0)*(R1-R0)*(1.0-e)+2.0/(R3+R2)*(R3-R2)*e)/pbstp/1e+5;
		if(maa<0) maa=0;
		/* Activation enthalpy recalc using enthalpy changes */
		/* Current Enthalpy */
		mhh1=((H0*(1.0-n)+H1*n)*(1.0-e)+(H2*(1.0-n)+H3*n)*e);
		/* Pmin Enthalpy */
		mhh0=(td[n1][0 ][mm3][1]*(1.0-e) + td[n1+1][0 ][mm3][1]*e)*1000.0*4.1837;
		/* Enthalpy Difference calc */
		mdhh=(mhh1-mhh0);

		/* Save TD variables */
		*Mgg=mgg;
		*Mro=mro;
		*Mwa=mwa;
		*Mcp=mcp;
		*Mbb=mbb;
		*Maa=maa;
		*Mdhh=mdhh;
		*Mkt+=krad;
	}

	/* Deep TD base type */
	if(1==0 || mpb>0.75*pbmax || ynpb==1)
	{
		switch (mm2)
		{
			/* MORB DATABASE */
			/* UPPER, LOWER Crust */
			case 5:
			case 6:
			case 17:
			case 18:
			case 37:
			case 38:
			/* Sediments */
			case 2:
			case 3:
			case 4:
			/* Molten Sediments */
			case 22:
			case 23:
			case 24:
			/* Molten crust */
			case 25:
			case 26:
			/* Basalt */
			case 16:
			case 7:
			/* Molten Basalt */
			case 36:
			case 27:
			/* Gabbro */
			case 8:
			/* Molten Gabbro */
			case 28: mm3=10; break;
			/**/
			/* PIROLITE DATABASE */
			/* Dry peridotite */
			case 9:
			case 12:
			case 14:
			case 10:
			/* Wet peridotite */
			case 13:
			case 11:
			/* Molten peridotite */
			case 34: mm3=9; break;
			// Added missing rock types
			case 15:
			case 19:
			case 20:
			case 21:
			case 29:
			case 30:
			/* Unknown type */
			default: {fprintf(fp_log,"Deep TD: Unknown rock type for TD database %d, for marker %ld with T= %f, P=%f \n",mm2,mm1,mtk,mpb); fflush(fp_log); exit(0);}
		}
		/* ABCD-4Cell Number */
		e=(mtk-tkmin1)/tkstp1;
		if(e<0) e=0;
		if(e>(double)(tknum1-1)) e=(double)(tknum1-1);
		n=(mpb-pbmin1)/pbstp1;
		if(n<0) n=0;
		if(n>(double)(pbnum1-1)) n=(double)(pbnum1-1);
		n1=(int)(e);
		if(n1>tknum1-2) n1=tknum1-2;
		n2=(int)(n);
		if(n2>pbnum1-2) n2=pbnum1-2;
		/* e,n Calc */
		e=(e-(double)(n1));
		n=(n-(double)(n2));
		/* Ro H values */
		/* 0 2 */
		/* 1 3 */
		R0=td[n1  ][n2  ][mm3][0]*1000.0;
		R1=td[n1  ][n2+1][mm3][0]*1000.0;
		R2=td[n1+1][n2  ][mm3][0]*1000.0;
		R3=td[n1+1][n2+1][mm3][0]*1000.0;
		H0=td[n1  ][n2  ][mm3][1]*1000.0*4.1837;
		H1=td[n1  ][n2+1][mm3][1]*1000.0*4.1837;
		H2=td[n1+1][n2  ][mm3][1]*1000.0*4.1837;
		H3=td[n1+1][n2+1][mm3][1]*1000.0*4.1837;
		W0=td[n1  ][n2  ][mm3][4];
		W1=td[n1  ][n2+1][mm3][4];
		W2=td[n1+1][n2  ][mm3][4];
		W3=td[n1+1][n2+1][mm3][4];
		G0=td[n1  ][n2  ][mm3][3]*1000.0;G0*=G0*R0;
		G1=td[n1  ][n2+1][mm3][3]*1000.0;G1*=G1*R1;
		G2=td[n1+1][n2  ][mm3][3]*1000.0;G2*=G2*R2;
		G3=td[n1+1][n2+1][mm3][3]*1000.0;G3*=G3*R3;
		/* Shear modulus calc by interpolation */
		mgg=((G0*(1.0-n)+G1*n)*(1.0-e)+(G2*(1.0-n)+G3*n)*e);
		/* Ro calc by interpolation */
		mro=((R0*(1.0-n)+R1*n)*(1.0-e)+(R2*(1.0-n)+R3*n)*e);
		/* Water wt% calc by interpolation */
		mwa=0;
		/* Water in crystals */
		if(mm2!=9 && mm2!=10 && mm2!=14 && mpb<235000.0) 
		{
			dmwa=0.1;
			mwa+=dmwa;
			wro=1050.0;
			mro=100.0/((100.0-dmwa)/mro+dmwa/wro);
		}
		/* Cp calc by interpolation */
		mcp=((H2-H0)*(1.0-n)+(H3-H1)*n)/tkstp1;
		if(mcp<1e+2) mcp=1e+2; else if(mcp>5e+4) mcp=5e+4;
		/* Effective adiabatic betta=1/V*dV/dT=ro/T*[-dH/dP+V] calc by interpolation */
		mbb=(2.0/(R1+R0)-(H1-H0)/pbstp1/1e+5)*(1.0-e)+(2.0/(R3+R2)-(H3-H2)/pbstp1/1e+5)*e;
		mbb*=mro/mtk;
		if(mbb<-1e-2) mbb=-1e-2; else if(mbb>1e-2) mbb=1e-2;
		/* Effective compressibility term alpha=1/ro*d(ro)/dP calc by interpolation */
		maa=(2.0/(R1+R0)*(R1-R0)*(1.0-e)+2.0/(R3+R2)*(R3-R2)*e)/pbstp1/1e+5;
		if(maa<0) maa=0;
		/* Activation enthalpy recalc using enthalpy changes */
		/* Current Enthalpy */
		mhh1=((H0*(1.0-n)+H1*n)*(1.0-e)+(H2*(1.0-n)+H3*n)*e);
		/* Pmin Enthalpy */
		mhh0=(td[n1][0 ][mm3][1]*(1.0-e) + td[n1+1][0 ][mm3][1]*e)*1000.0*4.1837;
		/* Enthalpy Difference calc */
		mdhh=(mhh1-mhh0);
		/* Thermal conductivity */
		mkt=mkt1+krad;

		/* Computing transitional parameters */
		if(1==0 || mpb>pbmax || ynpb==1)
			// Manny has 1==1
		{
			/* Save TD variables */
			*Mgg=mgg;
			*Mro=mro;
			*Mwa=mwa;
			*Mcp=mcp;
			*Mbb=mbb;
			*Maa=maa;
			*Mdhh=mdhh;
			*Mkt=mkt;
		}
		else
		{
			xold=(pbmax-mpb)/(0.25*pbmax);
			/* Save TD variables */
			// Second column comes from shallow database assignment above, but I never reach into this deep one ! 
			mgg=mgg*(1.0-xold)+ *Mgg *xold;
			mro=mro*(1.0-xold)+ *Mro *xold;
			mwa=mwa*(1.0-xold)+ *Mwa *xold;
			mcp=mcp*(1.0-xold)+ *Mcp *xold;
			mbb=mbb*(1.0-xold)+ *Mbb *xold;
			maa=maa*(1.0-xold)+ *Maa *xold;
			mdhh=mdhh*(1.0-xold)+ *Mdhh *xold;
			mkt=mkt*(1.0-xold)+ *Mkt *xold;
			*Mgg=mgg;
			*Mro=mro;
			*Mwa=mwa;
			*Mcp=mcp;
			*Mbb=mbb;
			*Maa=maa;
			*Mdhh=mdhh;
			*Mkt=mkt;
		}
	}
}
/* End OMP Thermodynamic database use for ro, Cp */

// *** Interpolation routines using the following nodal locations *** 
/* Staggered Nodes num */
/*   [0]                [3]                [6]   */
/*  T0,xy0    Vy0     T3,xy3     Vy3             */
/*                                               */
/*   Vx0    P4,xx4,yy4  Vx3    P7,xx7,yy7        */
/*                                               */
/*   [1]                [4]                [7]   */
/*  T,xy1     Vy1     T4,xy4     Vy4             */
/*                                               */
/*   Vx1    P5,xx5,yy5  Vx4    P8,xx8,yy8        */
/*                                               */
/*   [2]                [5]                [8]   */
/*                                               */
/*                                               */


/* Weights for horizontal and vertical nodes calculation for marker interpolation */ 
void nodewt(long int m1min, long int m1max, long int m2min, long int m2max, double x, double y, int ynx, int yny)
	/* m1min,m1max, m2min,m2max - node X,Y number limits */
	/* x,y - current pont coordinates */
	/* ynx, yny - Type of shifts: No(0), Back(-1), Forw(1) */
{
	/* Eyy vertical position */
	long int m3;
	int nx,ny;

	/* Weigths in horizontal directions */
	/* Load distances to xn[] */
	if(ynx<0) 
	{
		for (m3=m1min;m3<=m1max;m3++)
		{
			xn[m3-m1min]=(gx[m3]+gx[m3-1])/2.0;
		}
	}
	if(ynx==0) 
	{
		for (m3=m1min;m3<=m1max;m3++)
		{
			xn[m3-m1min]=gx[m3];
		}
	}
	if(ynx>0) 
	{
		for (m3=m1min;m3<=m1max;m3++)
		{
			xn[m3-m1min]=(gx[m3]+gx[m3+1])/2.0;
		}
	}

	/* Calc maximal position in xn[] */
	nx=(int)(m1max-m1min);

	/* Calc coefficients for horizontal direction */
	fdweight(nx,0,x);
	/**/
	/* Reload horizontal coefficients to cn[] */
	for (m3=0;m3<=nx;m3++)
	{
		cn[m3][1]=cn[m3][0];
	}
	
	/* Weigths in vertical directions */
	/* Load distances to xn[] */
	if(yny<0) 
	{
		for (m3=m2min;m3<=m2max;m3++)
		{
			xn[m3-m2min]=(gy[m3]+gy[m3-1])/2.0;
		}
	}
	if(yny==0) 
	{
		for (m3=m2min;m3<=m2max;m3++)
		{
			xn[m3-m2min]=gy[m3];
		}
	}
	if(yny>0) 
	{
		for (m3=m2min;m3<=m2max;m3++)
		{
			xn[m3-m2min]=(gy[m3]+gy[m3+1])/2.0;
		}
	}

	/* Calc maximal position in xn[] */
	ny=(int)(m2max-m2min);

	/* Calc coefficients for horizontal direction */
	fdweight(ny,0,y);
}
/* End Weights for horizontal and vertical nodes calculation for marker interpolation */ 



/* Calculation of EE,VX,VY,ESP, and PR by Interpolation */
void allinteriomp(double x, double y, long int m10, long int m20, double *VX, double *VY, double *PR, double *ESP, double *EE)
	/* x,y - XY location of point for Vx,Vy calc */
{
	/* Counters */
	long int m1,m2,m3;
	/* en-NormalisedDistance */
	// Keep EXX and EXY local here, so only calculates EE
	double e,n,ival,xrat,EXX,EXY;
	/**/
	/**/
	/* Check X,Y */
	if(x<0) x=0; else if(x>xsize) x=xsize;
	if(y<0) y=0; else if(y>ysize) y=ysize;
	/**/
	/**/
	/**/
	/* Check weighting for interpolation */
	xrat=2.0/3.0;
	if(x<(gx[0]+gx[1])/2.0) xrat=1.0;
	if(x>(gx[xnumx-2]+gx[xnumx-1])/2.0) xrat=1.0;
	if(y<(gy[0]+gy[1])/2.0) xrat=1.0;
	if(y>(gy[ynumy-2]+gy[ynumy-1])/2.0) xrat=1.0;
	/**/
	/**/
	/**/
	// Store for more usage throughout subroutine
	m1=m10;
	m2=m20;
	/**/
	/**/
	/**/	
	/* EXY, ESP interpolation ------------------------ */
	// Clear buffer
	*ESP=0;
	/* Horizontal,Vertical limits for interpolation calc */
	if(m10<1) m10=1; if(m10>xnumx-3) m10=xnumx-3;
	if(m20<1) m20=1; if(m20>ynumy-3) m20=ynumy-3;
	/**/
	/* Calc normalized distances */
	// Note that the nodal distance is now fixed, while before could change it with intermod. If want that again see old scripts in dynwif/CleanOldRun..
	e=(x-gx[m10])/(gx[m10+1]-gx[m10]);
	n=(y-gy[m20])/(gy[m20+1]-gy[m20]);
	/* Vx interpolation ------------------------ */
	m3=m10*ynumy+m20;
	EXY=(1.0-e)*(1.0-n)*exy[m3]+(1.0-e)*n*exy[m3+1]+e*(1.0-n)*exy[m3+ynumy]+e*n*exy[m3+ynumy+1];
	*ESP=(1.0-e)*(1.0-n)*esp[m3]+(1.0-e)*n*esp[m3+1]+e*(1.0-n)*esp[m3+ynumy]+e*n*esp[m3+ynumy+1];
	/* End EXY, ESP interpolation ------------------------ */
	/**/
	/**/
	/**/
	/* Exx, P interpolation ------------------------ */	
	// Reset and clear buffer
	m10=m1;
	m20=m2;
	*EE=0;
	*PR=0;
	*VX=0; 
	*VY=0;
	/* Horizontal,Vertical limits for interpolation calc */
	if(x>(gx[m10]+gx[m10+1])/2.0) m10++;
	if(y>(gy[m20]+gy[m20+1])/2.0) m20++;	
	if(m10<1) m10=1; if(m10>xnumx-2) m10=xnumx-2;
	if(m20<1) m20=1; if(m20>ynumy-2) m20=ynumy-2;
	/* Calc normalized distances */
	e=(x-(gx[m10-1]+gx[m10])/2.0)/((gx[m10+1]-gx[m10-1])/2.0);
	n=(y-(gy[m20-1]+gy[m20])/2.0)/((gy[m20+1]-gy[m20-1])/2.0);
	/* Interpolation ------------------------ */
	m3=m10*ynumy+m20;
	EXX=(1.0-e)*(1.0-n)*exx[m3]+(1.0-e)*n*exx[m3+1]+e*(1.0-n)*exx[m3+ynumy]+e*n*exx[m3+ynumy+1];	
	// QUESTION TARAS why Interpolate pressure here, do already in interp or d? Now I do port it back, so rm if no need ...	
	// I guess you could also formulate this more in general, as sometimes I have the feeling some variables are interpolated needlessly. Could you please go over these routines and removed what is not really need to speed the code up?
	*PR=(1.0-e)*(1.0-n)*pr[m3]+(1.0-e)*n*pr[m3+1]+e*(1.0-n)*pr[m3+ynumy]+e*n*pr[m3+ynumy+1];
	// Include small weight (xrat) from farther away nodes for velocities	
	*VX=( (1.0-e)*(1.0-n)*(vx[m3-1]+vx[m3-ynumy-1])+(1.0-e)*n*(vx[m3]+vx[m3-ynumy])             +e*(1.0-n)*(vx[m3+ynumy-1]+vx[m3-1])+e*n*(vx[m3+ynumy]+vx[m3]) ) * 0.5*(1.0-xrat);
	*VY=( (1.0-e)*(1.0-n)*(vy[m3-ynumy]+vy[m3-ynumy-1])+(1.0-e)*n*(vy[m3-ynumy+1]+vy[m3-ynumy]) +e*(1.0-n)*(vy[m3]+vy[m3-1])+e*n*(vy[m3+1]+vy[m3]) ) * 0.5*(1.0-xrat);	
	//eps[11]+=ival*(vx[m3-1]+vx[m3-ynumy-1])*0.5*(1.0-xrat); QUESTION TARAS Why use nodes above and left above here ??
	//eps[12]+=ival*(vy[m3-ynumy]+vy[m3-ynumy-1])*0.5*(1.0-xrat);
	// Calculate second invariant
	*EE=pow(EXX*EXX+EXY*EXY,0.5);	
	/* End SIGxx*EPSxx,SIGyy*EPSyy interpolation ------------------------ */

	/* Vx interpolation ------------------------ */
	// Reset and clear buffer
	m10=m1;
	m20=m2;
	/* Horizontal,Vertical limits for interpolation calc */
	if(y<(gy[m20]+gy[m20+1])/2.0) m20-=1;
	if(m10<0) m10=0; if(m10>xnumx-2) m10=xnumx-2;
	if(m20<0) m20=0; if(m20>ynumy-3) m20=ynumy-3;
	/* Calc normalized distances */
	e=(x-gx[m10])/(gx[m10+1]-gx[m10]);
	n=(y-(gy[m20]+gy[m20+1])/2.0)/((gy[m20+2]-gy[m20])/2.0);
	/* Vx interpolation ------------------------ */
	m3=m10*ynumy+m20;
	*VX+=((1.0-e)*(1.0-n)*vx[m3]+(1.0-e)*n*vx[m3+1]+e*(1.0-n)*vx[m3+ynumy]+e*n*vx[m3+ynumy+1])*xrat;
	/* End Vx interpolation ------------------------ */
	/**/
	/**/
	/**/
	/* Vy interpolation ------------------------ */
	// Reset and clear buffer
	m10=m1;
	m20=m2;
	/* Horizontal,Vertical limits for interpolation calc */
	if(x<(gx[m10]+gx[m10+1])/2.0) m10-=1;	
	if(m10<0) m10=0; if(m10>xnumx-3) m10=xnumx-3;
	if(m20<0) m20=0; if(m20>ynumy-2) m20=ynumy-2;
	/* Calc normalized distances */
	e=(x-(gx[m10]+gx[m10+1])/2.0)/((gx[m10+2]-gx[m10])/2.0);
	n=(y-gy[m20])/(gy[m20+1]-gy[m20]);
	/* Vy interpolation ------------------------ */
	m3=m10*ynumy+m20;
	*VY+=((1.0-e)*(1.0-n)*vy[m3]+(1.0-e)*n*vy[m3+1]+e*(1.0-n)*vy[m3+ynumy]+e*n*vy[m3+ynumy+1])*xrat;
	/* End Vy interpolation ------------------------ */
	/**/
	/**/
	/**/
}
/* OMP Interpolate Vx,Vy, EPSxx,EPSyy,EPSxy, SPINxy from surrounding nodes to marker at x,y  */



/* OMP Calculation of T,T0 for current location by Interpolation */
void allintertomp(double x, double y, long int m10, long int m20, double *TK, double *TK2)
	/* x,y - XY location of point for Vx,Vy calc */
	/* m10, m20 - Upper left node */
	// TK - marker temperature
{
	/* Counters */
	long int m3;
	/* en-NormalizedDistance */
	double e,n,ival;

	/* Check X,Y */
	if(x<0) x=0; else if(x>xsize) x=xsize;
	if(y<0) y=0; else if(y>ysize) y=ysize;

	/* T interpolation ------------------------ */
	/* Buffer clear */
	*TK=*TK2=0;
	
	/* Horizontal,Vertical limits for interpolation calc */
	if(m10<0) m10=0; if(m10>xnumx-2) m10=xnumx-2;
	if(m20<0) m20=0; if(m20>ynumy-2) m20=ynumy-2;
	
	/* Calc normalized distances */
	e=(x-gx[m10])/(gx[m10+1]-gx[m10]);
	n=(y-gy[m20])/(gy[m20+1]-gy[m20]);
	
	/* T interpolation ------------------------ */
	m3=m10*ynumy+m20;
	*TK=(1.0-e)*(1.0-n)*tk[m3]+(1.0-e)*n*tk[m3+1]+e*(1.0-n)*tk[m3+ynumy]+e*n*tk[m3+ynumy+1];
	*TK2=(1.0-e)*(1.0-n)*tk2[m3]+(1.0-e)*n*tk2[m3+1]+e*(1.0-n)*tk2[m3+ynumy]+e*n*tk2[m3+ynumy+1];
	
	/* End T interpolation ------------------------ */
}
/* OMP Calculation of T,T0 for current location by Interpolation */




/* OMP Calculation of P by Interpolation */
double allinterpomp(double x, double y, long int m10, long int m20)
	/* x,y - XY location of point for Vx,Vy calc */
	/* m10, m20 - Upper left node */
{
	/* Counters */
	long int m3;
	/* en-Normalized distance */
	double ival,e,n;
        
	/* Check X,Y */
	if(x<0) x=0; else if(x>xsize) x=xsize;
	if(y<0) y=0; else if(y>ysize) y=ysize;

	/* Buffer clear */
	ival=0;
	/* Horizontal,Vertical limits for interpolation calc */
	if(x>(gx[m10]+gx[m10+1])/2.0) m10++;
	if(y>(gy[m20]+gy[m20+1])/2.0) m20++;	
	if(m10<1) m10=1; if(m10>xnumx-2) m10=xnumx-2;
	if(m20<1) m20=1; if(m20>ynumy-2) m20=ynumy-2;
	/* Calc normalized distances */
	e=(x-(gx[m10-1]+gx[m10])/2.0)/((gx[m10+1]-gx[m10-1])/2.0);
	n=(y-(gy[m20-1]+gy[m20])/2.0)/((gy[m20+1]-gy[m20-1])/2.0);
	/* P interpolation ------------------------ */
	m3=m10*ynumy+m20;
	ival=(1.0-e)*(1.0-n)*pr[m3]+(1.0-e)*n*pr[m3+1]+e*(1.0-n)*pr[m3+ynumy]+e*n*pr[m3+ynumy+1];	
	/* Return pressure */
	return ival;	
	/*
	fprintf(fp_log,"eps %e %e ",m1,m2,e,n); getchar();
	if(timestep){fprintf(fp_log,"P1 %e %e  %ld %ld %e %e %e",x,y,m10,m20,e,n,ival);getchar();}
	*/
}
/* OMP End calculation of P by Interpolation */



/* Calculation of SIGij by Interpolation */
void allinterdomp(double x, double y,long int m10, long int m20, double *TK,double *EXY,double *EXYE,double *SXY,double *SXYE,double *EXX,double *SXX,double *PR,double *SXXE,double *SPPE,double *EXXE,double *VX, double *MVX, double *VY, double *MVY)
	/* x,y - XY location of point for Vx,Vy calc */
{
	/* Counters */
	long int m1,m2,m3;
	/* en-NormalisedDistance */
	double ival,e,n,xrat;
	/**/
	/**/
	/* Check X,Y */
	if(x<0) x=0; else if(x>xsize) x=xsize;
	if(y<0) y=0; else if(y>ysize) y=ysize;
	/**/
	/**/
	/**/
	/* Store Up Left Node X,Y Num for later re-usage */
	m1=m10;
	m2=m20;
	/**/
	/**/
	/* Check weighting for interpolation */
	xrat=2.0/3.0;
	if(x<(gx[0]+gx[1])/2.0) xrat=1.0;
	if(x>(gx[xnumx-2]+gx[xnumx-1])/2.0) xrat=1.0;
	if(y<(gy[0]+gy[1])/2.0) xrat=1.0;
	if(y>(gy[ynumy-2]+gy[ynumy-1])/2.0) xrat=1.0;
	/**/
	/**/
	/* T interpolation ------------------------ */
	/* Buffer clear */
	*TK=0;
	/* Horizontal,Vertical limits for interpolation calc */
	if(m10<0) m10=0; if(m10>xnumx-2) m10=xnumx-2;
	if(m20<0) m20=0; if(m20>ynumy-2) m20=ynumy-2;
	/* Calc normalized distances */
	e=(x-gx[m10])/(gx[m10+1]-gx[m10]);
	n=(y-gy[m20])/(gy[m20+1]-gy[m20]);
	/* EPSxy Interpolate after interpolation weights */
	m3=m10*ynumy+m20;
	*TK=(1.0-e)*(1.0-n)*tk[m3]+(1.0-e)*n*tk[m3+1]+e*(1.0-n)*tk[m3+ynumy]+e*n*tk[m3+ynumy+1];
	/**/
	/* End SIGij old interpolation ------------------------ */
	/* End T interpolation ------------------------ */
	/**/
	/**/
	/**/
	/* SIGxy interpolation ------------------------ */
	// Reset and clear buffer 
	m10=m1;
	m20=m2;
	*EXY=*EXYE=*SXY=*SXYE=0;
	/* Horizontal,Vertical limits for interpolation calc */
	if(m10<1) m10=1; if(m10>xnumx-3) m10=xnumx-3;
	if(m20<1) m20=1; if(m20>ynumy-3) m20=ynumy-3;
	/* Calc normalized distances */
	e=(x-gx[m10])/(gx[m10+1]-gx[m10]);
	n=(y-gy[m20])/(gy[m20+1]-gy[m20]);
	/**/
	/* EPSxy Interpolate after interpolation weights */
	m3=m10*ynumy+m20;
	*EXY=(1.0-e)*(1.0-n)*exy[m3]+(1.0-e)*n*exy[m3+1]+e*(1.0-n)*exy[m3+ynumy]+e*n*exy[m3+ynumy+1];
	*EXYE=(1.0-e)*(1.0-n)*exye[m3]+(1.0-e)*n*exye[m3+1]+e*(1.0-n)*exye[m3+ynumy]+e*n*exye[m3+ynumy+1];
	*SXY=(1.0-e)*(1.0-n)*sxy[m3]+(1.0-e)*n*sxy[m3+1]+e*(1.0-n)*sxy[m3+ynumy]+e*n*sxy[m3+ynumy+1];
	*SXYE=(1.0-e)*(1.0-n)*sxye[m3]+(1.0-e)*n*sxye[m3+1]+e*(1.0-n)*sxye[m3+ynumy]+e*n*sxye[m3+ynumy+1];
	/* End SIGxy interpolation ------------------------ */
	/**/
	/**/
	/**/
	/* SIGxx,SIGyy interpolation ------------------------ */
	// Reset and clear buffer 
	m10=m1;
	m20=m2;
	*EXX=*SXX=*PR=*SXXE=*SPPE=*EXXE=0;
	*VX=*MVX=0;
	*VY=*MVY=0;
	/* Horizontal,Vertical limits for interpolation calc */
	if(x>(gx[m10]+gx[m10+1])/2.0) m10++;
	if(y>(gy[m20]+gy[m20+1])/2.0) m20++;
	if(m10<1) m10=1; if(m10>xnumx-2) m10=xnumx-2;
	if(m20<1) m20=1; if(m20>ynumy-2) m20=ynumy-2;
	/* Calc normalized distances */
	e=(x-(gx[m10-1]+gx[m10])/2.0)/((gx[m10+1]-gx[m10-1])/2.0);
	n=(y-(gy[m20-1]+gy[m20])/2.0)/((gy[m20+1]-gy[m20-1])/2.0);
	/* P interpolation ------------------------ */
	m3=m10*ynumy+m20;
	*EXX =(1.0-e)*(1.0-n)*exx[m3]+(1.0-e)*n*exx[m3+1]+e*(1.0-n)*exx[m3+ynumy]+e*n*exx[m3+ynumy+1];
	*SXX =(1.0-e)*(1.0-n)*sxx[m3]+(1.0-e)*n*sxx[m3+1]+e*(1.0-n)*sxx[m3+ynumy]+e*n*sxx[m3+ynumy+1];
	*PR  =(1.0-e)*(1.0-n)*pr[m3]+(1.0-e)*n*pr[m3+1]+e*(1.0-n)*pr[m3+ynumy]+e*n*pr[m3+ynumy+1];
	*SPPE=(1.0-e)*(1.0-n)*sppe[m3]+(1.0-e)*n*sppe[m3+1]+e*(1.0-n)*sppe[m3+ynumy]+e*n*sppe[m3+ynumy+1];
	*SXXE=(1.0-e)*(1.0-n)*sxxe[m3]+(1.0-e)*n*sxxe[m3+1]+e*(1.0-n)*sxxe[m3+ynumy]+e*n*sxxe[m3+ynumy+1];
	*EXXE=(1.0-e)*(1.0-n)*exxe[m3]+(1.0-e)*n*exxe[m3+1]+e*(1.0-n)*exxe[m3+ynumy]+e*n*exxe[m3+ynumy+1];
	*VX=( (1.0-e)*(1.0-n)*(vx[m3-1]+vx[m3-ynumy-1])+(1.0-e)*n*(vx[m3]+vx[m3-ynumy])             +e*(1.0-n)*(vx[m3+ynumy-1]+vx[m3-1])+e*n*(vx[m3+ynumy]+vx[m3]) )*0.5 *(1.0-xrat);
	*VY=( (1.0-e)*(1.0-n)*(vy[m3-ynumy]+vy[m3-ynumy-1])+(1.0-e)*n*(vy[m3-ynumy+1]+vy[m3-ynumy]) +e*(1.0-n)*(vy[m3]+vy[m3-1])+e*n*(vy[m3+1]+vy[m3]) )*0.5 *(1.0-xrat);
	*MVX=( (1.0-e)*(1.0-n)*(mvx[m3-1]+mvx[m3-ynumy-1])+(1.0-e)*n*(mvx[m3]+mvx[m3-ynumy])             +e*(1.0-n)*(mvx[m3+ynumy-1]+mvx[m3-1])+e*n*(mvx[m3+ynumy]+mvx[m3]) )*0.5 *(1.0-xrat);
	*MVY=( (1.0-e)*(1.0-n)*(mvy[m3-ynumy]+mvy[m3-ynumy-1])+(1.0-e)*n*(mvy[m3-ynumy+1]+mvy[m3-ynumy]) +e*(1.0-n)*(mvy[m3]+mvy[m3-1])+e*n*(mvy[m3+1]+mvy[m3]) )*0.5 *(1.0-xrat);
	/* End SIGxx,SIGyy interpolation ------------------------ */
	/**/
	/**/
	/**/
	/* Vx interpolation ------------------------ */
	// Reset and clear buffer
	m10=m1;
	m20=m2;
	/* Horizontal,Vertical limits for interpolation calc */
	if(y<(gy[m20]+gy[m20+1])/2.0) m20-=1;
	if(m10<0) m10=0; if(m10>xnumx-2) m10=xnumx-2;
	if(m20<0) m20=0; if(m20>ynumy-3) m20=ynumy-3;
	/* Calc normalized distances */
	e=(x-gx[m10])/(gx[m10+1]-gx[m10]);
	n=(y-(gy[m20]+gy[m20+1])/2.0)/((gy[m20+2]-gy[m20])/2.0);
	/* Vx interpolation ------------------------ */
	m3=m10*ynumy+m20;
	//*VX=(1.0-e)*(1.0-n)*vx[m3]+(1.0-e)*n*vx[m3+1]+e*(1.0-n)*vx[m3+ynumy]+e*n*vx[m3+ynumy+1];
	//*MVX=(1.0-e)*(1.0-n)*mvx[m3]+(1.0-e)*n*mvx[m3+1]+e*(1.0-n)*mvx[m3+ynumy]+e*n*mvx[m3+ynumy+1];
	// Include small weight (xrat) from farther away nodes for velocities   
	*VX+=((1.0-e)*(1.0-n)*vx[m3]+(1.0-e)*n*vx[m3+1]+e*(1.0-n)*vx[m3+ynumy]+e*n*vx[m3+ynumy+1]) *xrat;
	*MVX+=((1.0-e)*(1.0-n)*mvx[m3]+(1.0-e)*n*mvx[m3+1]+e*(1.0-n)*mvx[m3+ynumy]+e*n*mvx[m3+ynumy+1]) *xrat;
	/* End Vx interpolation ------------------------ */
	/**/
	/**/
	/**/
	/* Vy interpolation ------------------------ */
	// Reset and clear buffer
	m10=m1;
	m20=m2;
	/* Horizontal,Vertical limits for interpolation calc */
	if(x<(gx[m10]+gx[m10+1])/2.0) m10-=1;	
	if(m10<0) m10=0; if(m10>xnumx-3) m10=xnumx-3;
	if(m20<0) m20=0; if(m20>ynumy-2) m20=ynumy-2;
	/* Calc normalized distances */
	e=(x-(gx[m10]+gx[m10+1])/2.0)/((gx[m10+2]-gx[m10])/2.0);
	n=(y-gy[m20])/(gy[m20+1]-gy[m20]);
	/* Vy interpolation ------------------------ */
	m3=m10*ynumy+m20;
	//*VY=(1.0-e)*(1.0-n)*vy[m3]+(1.0-e)*n*vy[m3+1]+e*(1.0-n)*vy[m3+ynumy]+e*n*vy[m3+ynumy+1];
	//*MVY=(1.0-e)*(1.0-n)*mvy[m3]+(1.0-e)*n*mvy[m3+1]+e*(1.0-n)*mvy[m3+ynumy]+e*n*mvy[m3+ynumy+1];
	// Include small weight (xrat) from farther away nodes for velocities   
	*VY+=((1.0-e)*(1.0-n)*vy[m3]+(1.0-e)*n*vy[m3+1]+e*(1.0-n)*vy[m3+ynumy]+e*n*vy[m3+ynumy+1]) *xrat;
	*MVY+=((1.0-e)*(1.0-n)*mvy[m3]+(1.0-e)*n*mvy[m3+1]+e*(1.0-n)*mvy[m3+ynumy]+e*n*mvy[m3+ynumy+1]) *xrat;
	/* End Vy interpolation ------------------------ */
	/*
	fprintf(fp_log,"eps %e %e ",m1,m2,e,n); getchar();
	*/
}
/* Calculation of SIGij by Interpolation */



/* Calculation of Vx,Vy, EPSxx*SIGxx,EPSyy*SIGyy,EPSxy*SIGxy by Interpolation */
// Not adapted for parallelization
void allinters(double x, double y)
	/* x,y - XY location of point for Vx,Vy calc */
{
	/* Counters */
	long int m1,m2,m3,m10,m20,m1min,m1max,m2min,m2max;
	/* en-NormalisedDistance */
	double ival;
	/**/
	/**/
	/* Check X,Y */
	if(x<0) x=0; else if(x>xsize) x=xsize;
	if(y<0) y=0; else if(y>ysize) y=ysize;
	/**/
	/**/
	/**/
	/* Up Left Node X,Y Num */
	wn[0]=m10=m1serch(x);
	wn[1]=m20=m2serch(y);
	/**/
	/**/
	/**/
	/* SIGxy*EPSxy interpolation ------------------------ */
	/* Buffer clear */
	eps[13]=0;
	/* Horizontal,Vertical limits for interpolation calc */
	m1min=m10; if(m1min<1) m1min=1; if(m1min>xnumx-3) m1min=xnumx-3;
	m1max=m1min+1+intermod; if(m1max>xnumx-2) m1max=xnumx-2;
	m1min=m1min-intermod; if(m1min<1) m1min=1;
	/**/
	m2min=m20; if(m2min<1) m2min=1; if(m2min>ynumy-3) m2min=ynumy-3;
	m2max=m2min+1+intermod; if(m2max>ynumy-2) m2max=ynumy-2;
	m2min=m2min-intermod; if(m2min<1) m2min=1;
	/**/
	/* Interpolation weights calc  after Fornberg (1996) */
	nodewt(m1min,m1max,m2min,m2max,x,y,0,0);
	/**/
	/* SIGxy,EPSxy Interpolate after interpolation weights */
	for (m1=m1min;m1<=m1max;m1++)
		for (m2=m2min;m2<=m2max;m2++)
	{
		/* Current node num, wt */
		m3=m1*ynumy+m2;
		ival=cn[m1-m1min][1]*cn[m2-m2min][0];
		eps[13]+=ival*sxy[m3]*sxy[m3]/(2.0*nu[m3]);
	}
	/* End SIGxy*EPSxy interpolation ------------------------ */
	/**/
	/**/
	/**/
	/* SIGxx*EPSxx, SIGyy*EPSyy interpolation ------------------------ */
	/* Buffer clear */
	eps[14]=0;
	/* Horizontal,Vertical limits for interpolation calc */
	m1min=m10; if(x>(gx[m10]+gx[m10+1])/2.0) m1min+=1;
	if(m1min<1) m1min=1; if(m1min>xnumx-2) m1min=xnumx-2;
	m1max=m1min+1+intermod; if(m1max>xnumx-1) m1max=xnumx-1;
	m1min=m1min-intermod; if(m1min<1) m1min=1;
	/**/
	m2min=m20; if(y>(gy[m20]+gy[m20+1])/2.0) m2min+=1;
	if(m2min<1) m2min=1; if(m2min>ynumy-2) m2min=ynumy-2;
	m2max=m2min+1+intermod; if(m2max>ynumy-1) m2max=ynumy-1;
	m2min=m2min-intermod; if(m2min<1) m2min=1;
	/**/
	/* Interpolation weights calc  after Fornberg (1996) */
	nodewt(m1min,m1max,m2min,m2max,x,y,-1,-1);
	/**/
	/* SIGxx,EPSxx,SIGyy,EPSyy,P Interpolate after interpolation weights */
	for (m1=m1min;m1<=m1max;m1++)
		for (m2=m2min;m2<=m2max;m2++)
	{
		/* Current node num, wt */
		m3=m1*ynumy+m2;
		ival=cn[m1-m1min][1]*cn[m2-m2min][0];
		eps[14]+=ival*sxx[m3]*sxx[m3]/(2.0*nd[m3]);
	}
	/* End SIGxx*EPSxx,SIGyy*EPSyy interpolation ------------------------ */
	/**/
	/**/
	/**/
	/* Vx interpolation ------------------------ */
	/* Buffer clear */
	eps[11]=0;
	/* Horizontal,Vertical limits for interpolation calc */
	m1min=m10-intermod; if(m1min<0) m1min=0;
	m1max=m10+1+intermod; if(m1max>xnumx-1) m1max=xnumx-1;
	/**/
	m2min=m20; if(y<(gy[m20]+gy[m20+1])/2.0) m2min-=1;
	if(m2min<0) m2min=0; if(m2min>ynumy-3) m2min=ynumy-3;
	m2max=m2min+1+intermod; if(m2max>ynumy-2) m2max=ynumy-2;
	m2min=m2min-intermod; if(m2min<0) m2min=0;
	/**/
	/* Interpolation weights calc  after Fornberg (1996) */
	nodewt(m1min,m1max,m2min,m2max,x,y,0,+1);
	/**/
	/* Vx Interpolate after interpolation weights */
	for (m1=m1min;m1<=m1max;m1++)
		for (m2=m2min;m2<=m2max;m2++)
	{
		/* Current node num, wt */
		m3=m1*ynumy+m2;
		ival=cn[m1-m1min][1]*cn[m2-m2min][0];
		eps[11]+=ival*vx[m3];
	}
	/* End Vx interpolation ------------------------ */
	/**/
	/**/
	/**/
	/* Vy interpolation ------------------------ */
	/* Buffer clear */
	eps[12]=0;
	/* Horizontal,Vertical limits for interpolation calc */
	m1min=m10; if(x<(gx[m10]+gx[m10+1])/2.0) m1min-=1;
	if(m1min<0) m1min=0; if(m1min>xnumx-3) m1min=xnumx-3;
	m1max=m1min+1+intermod; if(m1max>xnumx-2) m1max=xnumx-2;
	m1min=m1min-intermod; if(m1min<0) m1min=0;
	/**/
	m2min=m20-intermod; if(m2min<0) m2min=0;
	m2max=m20+1+intermod; if(m2max>ynumy-1) m2max=ynumy-1;
	/**/
	/* Interpolation weights calc  after Fornberg (1996) */
	nodewt(m1min,m1max,m2min,m2max,x,y,+1,0);
	/**/
	/* Vy Interpolate after interpolation weights */
	for (m1=m1min;m1<=m1max;m1++)
		for (m2=m2min;m2<=m2max;m2++)
	{
		/* Current node num, wt */
		m3=m1*ynumy+m2;
		ival=cn[m1-m1min][1]*cn[m2-m2min][0];
		eps[12]+=ival*vy[m3];
	}
	/* End Vy interpolation ------------------------ */
	/*
	fprintf(fp_log,"eps %e %e ",m1,m2,e,n); getchar();
	*/
}
/* Calculation of Vx,Vy, EPSxx*SIGxx,EPSyy*SIGyy,EPSxy*SIGxy by Interpolation */


/* Calculation of Vx,Vy, EPSxx,EPSyy,EPSxy, SPINxy by Interpolation */
// Not adapted for parallelization
void allinteri(double x, double y)
	/* x,y - XY location of point for Vx,Vy calc */
{
	/* Counters */
	long int m1,m2,m3,m10,m20,m1min,m1max,m2min,m2max;
	/* en-NormalisedDistance */
	double ival,xrat;
	/**/
	/**/
	/* Check X,Y */
	if(x<0) x=0; else if(x>xsize) x=xsize;
	if(y<0) y=0; else if(y>ysize) y=ysize;
	/**/
	/**/
	/**/
	/* Check weighting for interpolation */
	xrat=2.0/3.0;
	if(x<(gx[0]+gx[1])/2.0) xrat=1.0;
	if(x>(gx[xnumx-2]+gx[xnumx-1])/2.0) xrat=1.0;
	if(y<(gy[0]+gy[1])/2.0) xrat=1.0;
	if(y>(gy[ynumy-2]+gy[ynumy-1])/2.0) xrat=1.0;
	/**/
	/**/
	/**/
	/* Up Left Node X,Y Num */
	wn[0]=m10=m1serch(x);
	wn[1]=m20=m2serch(y);
	/**/
	/**/
	/**/
	/* EPSxy, SPINxy interpolation ------------------------ */
	/* Buffer clear */
	eps[4]=eps[30]=0;
	/* Horizontal,Vertical limits for interpolation calc */
	m1min=m10; if(m1min<1) m1min=1; if(m1min>xnumx-3) m1min=xnumx-3;
	m1max=m1min+1+intermod; if(m1max>xnumx-2) m1max=xnumx-2;
	m1min=m1min-intermod; if(m1min<1) m1min=1;
	/**/
	m2min=m20; if(m2min<1) m2min=1; if(m2min>ynumy-3) m2min=ynumy-3;
	m2max=m2min+1+intermod; if(m2max>ynumy-2) m2max=ynumy-2;
	m2min=m2min-intermod; if(m2min<1) m2min=1;
	/**/
	/* Interpolation weights calc  after Fornberg (1996) */
	nodewt(m1min,m1max,m2min,m2max,x,y,0,0);
	/**/
	/* SIGxy,EPSxy Interpolate after interpolation weights */
	for (m1=m1min;m1<=m1max;m1++)
		for (m2=m2min;m2<=m2max;m2++)
	{
		/* Current node num, wt */
		m3=m1*ynumy+m2;
		ival=cn[m1-m1min][1]*cn[m2-m2min][0];
		eps[4]+=ival*exy[m3];
		eps[30]+=ival*esp[m3];
	}
	/* End SIGxy*EPSxy interpolation ------------------------ */
	/**/
	/**/
	/**/
	/* EPSxx, EPSyy, P interpolation ------------------------ */
	/* Buffer clear */
	eps[6]=eps[10]=eps[11]=eps[12]=0;
	/* Horizontal,Vertical limits for interpolation calc */
	m1min=m10; if(x>(gx[m10]+gx[m10+1])/2.0) m1min+=1;
	if(m1min<1) m1min=1; if(m1min>xnumx-2) m1min=xnumx-2;
	m1max=m1min+1+intermod; if(m1max>xnumx-1) m1max=xnumx-1;
	m1min=m1min-intermod; if(m1min<1) m1min=1;
	/**/
	m2min=m20; if(y>(gy[m20]+gy[m20+1])/2.0) m2min+=1;
	if(m2min<1) m2min=1; if(m2min>ynumy-2) m2min=ynumy-2;
	m2max=m2min+1+intermod; if(m2max>ynumy-1) m2max=ynumy-1;
	m2min=m2min-intermod; if(m2min<1) m2min=1;
	/**/
	/* Interpolation weights calc  after Fornberg (1996) */
	nodewt(m1min,m1max,m2min,m2max,x,y,-1,-1);
	/**/
	/* SIGxx,EPSxx,SIGyy,EPSyy,P Interpolate after interpolation weights */
	for (m1=m1min;m1<=m1max;m1++)
		for (m2=m2min;m2<=m2max;m2++)
	{
		/* Current node num, wt */
		m3=m1*ynumy+m2;
		ival=cn[m1-m1min][1]*cn[m2-m2min][0];
		eps[6]+=ival*exx[m3];
		eps[10]+=ival*pr[m3];
		eps[11]+=ival*(vx[m3-1]+vx[m3-ynumy-1])*0.5*(1.0-xrat);
		eps[12]+=ival*(vy[m3-ynumy]+vy[m3-ynumy-1])*0.5*(1.0-xrat);
	}
	/* End SIGxx*EPSxx,SIGyy*EPSyy interpolation ------------------------ */
	/*
	depthp(x,y);eps[10]=eps[50];
	*/
	/**/
	/**/
	/**/
	/* Vx interpolation ------------------------ */
	/* Horizontal,Vertical limits for interpolation calc */
	m1min=m10-intermod; if(m1min<0) m1min=0;
	m1max=m10+1+intermod; if(m1max>xnumx-1) m1max=xnumx-1;
	/**/
	m2min=m20; if(y<(gy[m20]+gy[m20+1])/2.0) m2min-=1;
	if(m2min<0) m2min=0; if(m2min>ynumy-3) m2min=ynumy-3;
	m2max=m2min+1+intermod; if(m2max>ynumy-2) m2max=ynumy-2;
	m2min=m2min-intermod; if(m2min<0) m2min=0;
	/**/
	/* Interpolation weights calc  after Fornberg (1996) */
	nodewt(m1min,m1max,m2min,m2max,x,y,0,+1);
	/**/
	/* Vx Interpolate after interpolation weights */
	for (m1=m1min;m1<=m1max;m1++)
		for (m2=m2min;m2<=m2max;m2++)
	{
		/* Current node num, wt */
		m3=m1*ynumy+m2;
		ival=cn[m1-m1min][1]*cn[m2-m2min][0];
		eps[11]+=ival*vx[m3]*xrat;
	}
	/* End Vx interpolation ------------------------ */
	/**/
	/**/
	/**/
	/* Vy interpolation ------------------------ */
	/* Horizontal,Vertical limits for interpolation calc */
	m1min=m10; if(x<(gx[m10]+gx[m10+1])/2.0) m1min-=1;
	if(m1min<0) m1min=0; if(m1min>xnumx-3) m1min=xnumx-3;
	m1max=m1min+1+intermod; if(m1max>xnumx-2) m1max=xnumx-2;
	m1min=m1min-intermod; if(m1min<0) m1min=0;
	/**/
	m2min=m20-intermod; if(m2min<0) m2min=0;
	m2max=m20+1+intermod; if(m2max>ynumy-1) m2max=ynumy-1;
	/**/
	/* Interpolation weights calc  after Fornberg (1996) */
	nodewt(m1min,m1max,m2min,m2max,x,y,+1,0);
	/**/
	/* Vy Interpolate after interpolation weights */
	for (m1=m1min;m1<=m1max;m1++)
		for (m2=m2min;m2<=m2max;m2++)
	{
		/* Current node num, wt */
		m3=m1*ynumy+m2;
		ival=cn[m1-m1min][1]*cn[m2-m2min][0];
		eps[12]+=ival*vy[m3]*xrat;
	}
	/* End Vy interpolation ------------------------ */
	/*
	fprintf(fp_log,"eps %e %e ",m1,m2,e,n); getchar();
	*/
}
/* Calculation of Vx,Vy, EPSxx,EPSyy,EPSxy, SPINxy by Interpolation */


