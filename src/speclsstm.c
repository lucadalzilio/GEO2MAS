/* === Additional input parameters potentially needed on the fly === */
void flyinput_dynw()
{
	int mm1,im,ir;
	double gps_spacing;
	
	// --- Define additional or replacing input data ---
	// This will for once overwrite input read from prn (and set originally in init.t3c)
	
	// - Frictional parameters -
	// static_fric here replaces markb0/1[7] from init.t3c since need to be able to adapt through time

	// --- STM stage, w. events ---

	// Velocity Weakening
	//mus_mtvw = 0.5;			// Mix Di Toro et al. (2012), Den Hartog et al. (2012)
	mgamma_vw = 0.7;               		// = 1 - mu_d / mu_s , but also (beta-alpha)/mu_s in eq 3 in Ampuero & Ben-Zion (2008) -> mu_d = ( 1 - part_dyn ) * mu_s
	mvc_vw = 1.4e-9;                 	// Characteristic velocity, determining at which velocity friction is reduced by a half

	// Velocity Strengthening (only allowed on megathrust in PhD thesis van Dinther) 
	mus_vs = 0.35;
	mgamma_vs = -1.5;               	// = 1 - mu_d / mu_s , but also (beta-alpha)/mu_s in eq 3 in Ampuero & Ben-Zion (2008) -> mu_d = ( 1 - part_dyn ) * mu_s
	mvc_vs = 2e-9;              		// Characteristic velocity, determining at which velocity friction is reduced
	
	// lambfld = Pf/Ps = pore fluid pressure ratio (so not 1-this as in older version!)
	lambfld_stm = 0.0; //0.95;
			
	// Define location updip limit seismogenic zone in celsius
	tk_updipsez0 = 0 + 273;
	tk_updipsez1 = 0 + 273;   

	// --- Thermo-Mechanical stage, w.o. events ---
	mus_mtvw = 0.1;
	//mus_ini = 0.15;
	lambfld_ini = 0.999;
	
	// - Other parameters -
	// One-time shift HR grid to right by this many km, though keep this value to adjust hard-coded limits always	
	// PAY ATTENTION: to 0 when not shift done in simulation thusfar 
	shift_km = 0;
	
	// - Locations - 
	// nr of cells equal to 1 cm = 6.4 km above thrust interface as in lab comparison
	above_thrust = (int)(6400/(res_high));

	/* - - - pick events - - - */

	/* velocity threshold */
	VP_thresh = 5.0e-9;
	/* domain */
	X_11 = 1700e+3;
	X_21 = 2300e+3;
	Y_12 = 120e+3;

	// -------------------------- Select markers to be tracked -------------------------------------------
	
	// --- Surface (GPS) markers ---
	
	// GPS array (nr in head.c)
	gps_spacing = 10e3;		// km
	
	// Get locations for GPS surface markers assuming a line at t=0
	for (im=0;im<nr_gpsmarkers;im++)
	{
		xgps[im] = 1650e3+gps_spacing*im;
		// determine y from rock composition of not water or air from the top
	}
}

/* === Set global variables for stm depending on the time step === */
void set_stmsw()
{
	// 101 years and more,i.e. , 1000
	if (maxtmstep>100*3.15576e+7)
	{
		if (printmod) fprintf(fp_log,"\nLarge timestep part of the simulation... \n"); fflush(fp_log);

		inertyn = 0;				// Inertial term: yes (1) or no (0)
		start_cond = 0;				// Start post processing output: yes (1) or no (0)
		veldepfric = 0;			// Velocity Dependent Friction: yes (1) or no (0)

		// Lower strength at start when no slip velocity dependent weakening
                
		/* friction coefficients */
		/* basalt(7) and gabros(8) - oceanic crust */
		//markb0[7] = markb1[7] = 0.1;
		//markb0[8] = markb1[8] = 0.1;
	}
	// 7-100 years: increase strength interface
	if (maxtmstep>=7*3.15576e+7 && maxtmstep<=100*3.15576e+7)
	{
		if (printmod) fprintf(fp_log,"\nDecreasing timesteps... \n"); fflush(fp_log);

		inertyn = 0;
		start_cond = 0;
		veldepfric = 1;

		/* friction coefficients */
		/* basalt(7) and gabros(8) - oceanic crust */
		//markb0[7] = markb1[7] = 0.1;
		//markb0[8] = markb1[8] = 0.1;

		//lambfld = lambfld_stm;

		// To simplify setup plasticity only allowed at megathrust (=rock type 7) for first large-scale paper (van Dinther et al, 2013b, JGR)
#if setup==10
		// For all other normal rocktypes make both strength terms 0 (for markll through 1-1), so completely unabled:
		for (int ir=2;ir<7;ir++) {markll[ir] = 1; marka0[ir] = 0; marka1[ir] = 0;}
		for (int ir=8;ir<20;ir++) {markll[ir] = 1; marka0[ir] = 0; marka1[ir] = 0;}
		// Replaces 200 MPa cohesion, since was not always sufficient if changed setup / material parameters
#endif
	}
	// 6 years: transition
        if (maxtmstep>5*3.15576e+7 && maxtmstep<7*3.15576e+7)
        {
		if (printmod) fprintf(fp_log,"\nSeismic cycle timestepping transition... \n"); fflush(fp_log);	
	        inertyn = 1;
                start_cond = 0;
                veldepfric = 1;
	}
	// <=5 years: analyse runs
	if (maxtmstep<=5*3.15576e+7)
	{
		if (printmod) fprintf(fp_log,"\nSeismic cycle timestepping... \n"); fflush(fp_log);

		inertyn = 1;
		start_cond = 1;
		veldepfric = 1;

		// Keep high strength interface
		//lambfld = lambfld_stm;

		/* friction coefficients */
		/* basalt(7) and gabros(8) - oceanic crust */
		//markb0[7] = markb1[7] = 0.1;
		//markb0[8] = markb1[8] = 0.1;

		// To simplify setup plasticity only allowed at megathrust (=rock type 7) for first large-scale paper (van Dinther et al, 2013b, JGR)
#if setup==10
		// For all other normal rocktypes make both strength terms 0 (for markll through 1-1), so completely unabled:
		for (int ir=2;ir<7;ir++) {markll[ir] = 1; marka0[ir] = 0; marka1[ir] = 0;}
		for (int ir=8;ir<20;ir++) {markll[ir] = 1; marka0[ir] = 0; marka1[ir] = 0;}
		// Replaces 200 MPa cohesion, since was not always sufficient if changed setup / material parameters
#endif
	}
}


/* === Prepare post-processing analysis large-scale models === */
void preppostanalysis(double **Es1, int *nr_areas)
{
	long int m1;
	
	// Flag which markers need to be saved for picking algorithm, so keep fixed sequence in markers in hdf5
	for (m1=0;m1<marknum;m1++)
	{
		// Note is also hard-coded in mark.c !
		if ( marky[m1]<85e3 && markx[m1]>gx[m10_hr] && markx[m1]<gx[m11_hr] && markt[m1]>1 && markt[m1]<100)
		{
			// Rock marker (or molten)
			if ( markt[m1]<50)
			{
				follow[m1]=1;
				nm++;
			}
			// Fluid marker
			else 
			{
				follow[m1]=2;
				nmf++;
			}
		}
	}

	// Do not store follow array in prn, will reselect follow at dt=1
	// Better to reevaluate which markers to follow for evaluation of bending area markers that are close to HR limit
	// Though also need to post-process in two steps !! will see if necessary from time stamps hdf5 files...

	nm_old = nm;
	nmf_old = nmf;

	// Allocate memory dynamically once, can be resized if necessary by realloc
	mcf = (int*)malloc( sizeof(int) * nm * 1 );
	markxf = (double*)malloc( sizeof(double) * nm * 1 );
	markyf = (double*)malloc( sizeof(double) * nm * 1 );
	markef = (double*)malloc( sizeof(double) * nm * 1 );
	markexxf = (double*)malloc( sizeof(double) * nm * 1 );
	markexyf = (double*)malloc( sizeof(double) * nm * 1 );
	markxxf = (double*)malloc( sizeof(double) * nm * 1 );
	markxyf = (double*)malloc( sizeof(double) * nm * 1 );
	markvf = (double*)malloc( sizeof(double) * nm * 1 );
	markpf = (double*)malloc( sizeof(double) * nm * 1 );

	if (savefluid==1)
	{	
		mcff = (int*)malloc( sizeof(int) * nmf * 1 );
		markxff = (double*)malloc( sizeof(double) * nmf * 1 );
		markyff = (double*)malloc( sizeof(double) * nmf * 1 );
		markwff = (double*)malloc( sizeof(double) * nmf * 1 );
	}

	// Dynamically allocate memory to store total dissipated energy for all areas
	*nr_areas = 5;	
	*Es1 = (double*)malloc( sizeof(double) * *nr_areas * 1 );
}


/* === Post-processing of output at each timestep for printing to text files === */
void postproc_dynw(int *nr_areas, double *Es1, int f0)
{
	/* Counters */
	int m1,m2,m3,m4,irt,iy,na,count,pos,ii,*ncomp,*rocktype,index,pred_index,xresol_vis,yresol_vis,opt,*ythrust,*addythrust,im,mm1,imt;
	int diff_m1,diff_m2,m1_trench=0,m2_trench,m1_deeptrench,m2_deeptrench,m1_bottomwedge,m2_bottomwedge=1,deepest_bea,ok_rt,nrt,xend,ybegin,yend,localdiprange,xythrust;
	int m10_ca=0,m11_ca=0,m20_ca=0,m21_ca=0;
	double c10_ca,c11_ca,c20_ca,c21_ca,*maxstrain,*vpar_th,*vslip_th,*vslip_ne_th,*vx_th,*vy_th,*localdip,*gx_th,*gy_th,*sii_th,*nu_th,*gy_ath,*vpar_ath,*sii_ath,*nu_ath,*gy_bth,*vpar_bth,*sii_bth,*nu_bth;
	double area,EsDt,pra,siia,siia_nd,sxxa,sxya,eiia,exx_ne,exxa_ne,exy_ne,exya_ne,eii_ne,eiia_ne,shh_xx,shh_xy,dx,dy,timesum05,dt,xlim_phys,ylim_phys,xlim_DP_phys,ylim_DP_phys,xlim_UP_phys,ylim_UP_phys;
	char* numbers = ".";
	char fileDissEs[50];
	size_t szint,szfloat,szdouble;
	FILE *fld,*fldT,*fldT2,*fldT3,*flgps,*flphys,*DP_flphys,*UP_flphys,*flPK,*flCR;
	
	if (printmod) fprintf(fp_log,"Saving data each time step to *txt ...\n"); fflush(fp_log);
	
	count1++;
	rocktype = 0;
	// Set all elements of sii to zero, so can check if access the right elements
	for (ii=0; ii<MAXNOD; ii++) { sii[ii] = 0;}
	
	// - Calculate stress and strain rate second invariants -
	// Have deviatoric stresses sxx, sxy, to calculate visco-plastic strain rate            
	for (m1=m10_hr;m1<m11_hr+1;m1++)         // Each loop index end in C with one more, so also do the last index      
	{
		for (m2=m20_hr;m2<m21_hr+1;m2++)
		{   	// Get node number in 1-D array
			m3=m1*ynumy+m2;
			eii[m3]=sqrt(pow(exx[m3],2)+pow((exy[m3]+exy[m3+1]+exy[m3+ynumy]+exy[m3+ynumy+1])/4,2));
			sii[m3]=sqrt(pow(sxx[m3],2)+pow((sxy[m3]+sxy[m3+1]+sxy[m3+ynumy]+sxy[m3+ynumy+1])/4,2));
			sii_nd[m3]=sqrt(pow(sxy[m3],2)+pow(((sxx[m3]+pr[m3])+(sxx[m3+1]+pr[m3+1])+(sxx[m3+ynumy]+pr[m3+ynumy])+(sxx[m3+ynumy+1]+pr[m3+ynumy+1]))/4,2));
			sxx_nd[m3]= sxx[m3] - pr[m3];                                   //Are located all on the same node
			// sxy_nd = sxy, since off-diagonal components do not contain pressure component
		}
	}
	
	// - Automatically determine area(s) to perform calculations on -> limiting nodal points in arrays -
	// Interpolate rock types from markers to nodes; 0 = not compressed
	// index = nr visual nodes, determined in interpoler()
	xresol_vis = ((xsize/res_high)+1); 
	yresol_vis = ((ysize/res_high)+1); 
	pred_index = xresol_vis * yresol_vis;  
	ncomp = (int*)malloc( sizeof(int) * pred_index * 1 );
	opt = 5; 													//  5 means rock type, 60 is file number, but is not used in current version of interpoler()
	index = interpoler(ncomp, opt, 60, 0);
	// Correct for higher vis res in ncomp in low res part model
	diff_m1 = (int)(gx[m10_hr]/res_high) - m10_hr; 
	diff_m2=0;

	/* check position */
	//fprintf(fp_log,"++++++++++ m10_hr = %f, m11_hr = %f \n",gx[m10_hr]/1e3,gx[m11_hr]/1e3);
	//fprintf(fp_log,"++++++++++ m20_hr = %f, m21_hr = %f \n",gy[m20_hr]/1e3,gy[m21_hr]/1e3);
	
	/* - - - Get bottom wedge location - - - */
	for (m1=m10_hr;m1<m11_hr+1;m1++)
	{            
		for (m2=m20_hr;m2<m21_hr+1;m2++)
		{
			m3=(m1+diff_m1)*yresol_vis+(m2+diff_m2);
			if (ncomp[m3]==15) 
			{ 
				if(m2_bottomwedge<=m2)
				{
					m1_bottomwedge = m1; 	
					m2_bottomwedge = m2;
					m20_hr = m2;
					//fprintf(fp_log,"- - - m20_hr= %f m1_bottomwedge= %f, m2_bottomwedge= %f \n",gy[m20_hr]/1e3,gx[m1_bottomwedge]/1e3,gy[m2_bottomwedge]/1e3);
				}
				break;
			}
			if (ncomp[m3]!=0 && ncomp[m3]!=1 && ncomp[m3]!=3 && ncomp[m3]!=4 && ncomp[m3]!=5)
			{break;}
		}
	}
    /* reset value */
    m20_hr = 0;
    
	fprintf(fp_log,"- - - Find bottom wedge: is at x = %f, y = %f \n",gx[m1_bottomwedge]/1e3,gy[m2_bottomwedge]/1e3);
	check_bound(m10_hr,m11_hr,m1_bottomwedge,"postproc m1_bottomwedge");
	check_bound(m20_hr,m21_hr,m2_bottomwedge,"postproc m2_bottomwedge");
	fflush(fp_log);

	/* +++++++ Get trench location ++++++++ */
	/* define x direction --> check along y the rock type */

        for (m1=m10_hr;m1<m11_hr+1;m1++)
        {
	 	for (m2=m20_hr;m2<m21_hr+1;m2++)
		{
			m3=(m1+diff_m1)*yresol_vis+(m2+diff_m2);
			// Store and break if find air - or sediments? ---> try to check?
			// if sediments subduct at the trench, try to use rock type 3 (or 4?)
			if (ncomp[m3]==3) /* --> sediments at the trench? */
			{
				m1_trench=m1+20;
				m2_trench=m2-1;
				check_bound(m10_hr,m11_hr,m1_trench,"postproc m1_trench");
				check_bound(m20_hr,m21_hr,m2_trench,"postproc m2_trench");
				fprintf(fp_log,"Trench is at x = %f, y = %f \n",gx[m1_trench]/1e3,gy[m2_trench]/1e3);
				fflush(fp_log);
				break;
			}
		}
		// Also go out of this nested loop
		if (m1_trench!=0){break;}
	}

	timesum05 = timesum + 0.5*timestep;

	/* ---- START: DISSIPATED STRAIN ENERGY RATE ---- */
	for (na=0;na<*nr_areas;na++)
	{
		// Intialize variables for averages, for each area
		count = 0; EsDt = 0; pra = 0; siia = 0; siia_nd = 0; sxxa = 0; sxya = 0; eiia = 0; exxa_ne = 0; exya_ne = 0; eiia_ne = 0;
		
		// Determine area(s) to perform calculations on -> limiting nodal points in arrays
		// X: m10 ---- m11
		// Z: m20 |||| m21
		switch (na)
		{
			case 0:
				//Thrust              
				nrt = 2;
				rocktype    = (int*)malloc( sizeof(int) * nrt * 1 );
				rocktype[0] = 7;
				rocktype[1] = 8;
				// Would also like near by parts of rt 9, since this will more see transferance ...
				m10_ca = m1_trench;
				m11_ca = m1_trench+(int)(380e3/res_high);			// Wide, will be limited by rocktype
				m20_ca = m2_trench; 	       			  //(int)(9e3/res_high);		// Shallowest seismogenic zone of the world; 4 km below trench
				m21_ca = (int)(120e3/res_high);			  //(int)(90e3/res_high); 	// Deepest seismogenic zone of the world: 74km?  (Flores, Heuret et al. (2011))
				fprintf(fp_log,">>> Main-thrust area dimesions are: X= %f to %f, Y= %f to %f \n",gx[m10_ca],gx[m11_ca],gy[m20_ca],gy[m21_ca]);	
				break;
			case 1:
				// Bending area
				nrt = 4;
				rocktype    = (int*)malloc( sizeof(int) * nrt * 1 );
				rocktype[0] = 5;
				rocktype[1] = 6;
				rocktype[2] = 15;
				rocktype[3] = 16;
				m10_ca = m10_hr;
				m11_ca = m1_bottomwedge+(int)(50e3/res_high);
				m20_ca = (int)(15e3/res_high);
				m21_ca = (int)(100e3/res_high); 
				fprintf(fp_log,">>> Bending area dimesions are: X= %f to %f, Y= %f to %f \n",gx[m10_ca],gx[m11_ca],gy[m20_ca],gy[m21_ca]);	
				break;
			case 2:
				// Wedge
				nrt = 2;
				rocktype    = (int*)malloc( sizeof(int) * nrt * 1 );
				rocktype[0] = 3;
				rocktype[1] = 4;
				m10_ca = m1_trench;
				m11_ca = m1_bottomwedge+(int)(100e3/res_high);	// m1_trench+(int)(200e3/res_high) ;
				m20_ca = (int)(10e3/res_high);			// Wide, will be limited by rocktype
				m21_ca = m2_bottomwedge; 
				fprintf(fp_log,">>> Sedimentary wedge dimesions are: X= %f to %f, Y= %f to %f \n",gx[m10_ca],gx[m11_ca],gy[m20_ca],gy[m21_ca]);	
				break;
			case 3:
				// Overriding plate 
				nrt = 3;
				rocktype    = (int*)malloc( sizeof(int) * nrt * 1 );
				rocktype[0] = 13;
				rocktype[1] = 14;
				rocktype[2] = 17;
				rocktype[3] = 18;
				m10_ca = m1_bottomwedge;
				m11_ca = m1_bottomwedge+(int)(250e3/res_high);
				m20_ca = (int)(10e3/res_high);			// Wide, will be limited by rocktype
				m21_ca = (int)(100e3/res_high);			//(int)(25e3/res_high);	// Wide, will be limited by rocktype 
				fprintf(fp_log,">>> Upper plate area dimesions are: X= %f to %f, Y= %f to %f \n",gx[m10_ca],gx[m11_ca],gy[m20_ca],gy[m21_ca]);	
				break;
			// Also do it the old fashioned way, without rocktypes; +: less chaotic (as plasticity inherently is sensitive..)
			case 4:
				// All around thrust 
				// 100 means, no restrictions on rock type
				nrt    = 100;
				m10_ca = m10_hr;
				m11_ca = m1_bottomwedge+(int)(250e3/res_high);
				m20_ca = (int)(10e3/res_high);;
				m21_ca = (int)(120e3/res_high); 
				fprintf(fp_log,">>> All around orogenic belt: X= %f to %f, Y= %f to %f \n",gx[m10_ca],gx[m11_ca],gy[m20_ca],gy[m21_ca]);	
				fflush(fp_log);
				break;
		}
		
		// Calculate strain energy rate for each node, and then integrate over volume
		// In c also have space at the i=0 location !
		for (m1=m10_ca;m1<m11_ca+1;m1++)         // Each loop index end in C with one more, so also do the last index      
		{
			for (m2=m20_ca;m2<m21_ca+1;m2++)
			{       
				m4=(m1+diff_m1)*yresol_vis+(m2+diff_m2);
				
				// Check for array of rocktypes 
				ok_rt=0;
				if (nrt==100)
					{ok_rt=1; }
				else
				{
					for(irt=0;irt<nrt;irt++)
					{
						if (ncomp[m4]==rocktype[irt])
							{ok_rt = 1;	break;}
					}
				}
						
				// Only for certain rocktypes
				if (ok_rt==1)
				{
					count++;
					m3=m1*ynumy+m2;
					
					// Calculate viscous (and plastic?) shear heating shh for each node: strain rate e * stress s
					// Use normal, Non-Deviatoric stresses
					exx_ne = (sxx_nd[m3]/(2*nd[m3]));
					exy_ne = (sxy[m3]/(2*nu[m3]));
					eii_ne = sqrt(pow(exy_ne,2) + pow(exx_ne,2));
					shh_xx = exx_ne * sxx_nd[m3];
					shh_xy = exy_ne * sxy[m3];

					// From here add vals to right area 
					// Calculate Average stress values
					pra = pra + pr[m3];
					siia = siia + sii[m3];
					siia_nd = siia_nd + sii_nd[m3];
					sxxa = sxxa + sxx[m3];
					sxya = sxya + sxy[m3];
					
					// Calculate average non-elastic (ne) strain rates
					eiia    = eiia + eii[m3];
					eiia_ne = eiia_ne + eii_ne;
					exxa_ne = exxa_ne + exx_ne;
					exya_ne = exya_ne + exy_ne;
					
					// Integrate over volume V ( = dx*dy*1 ) -> Total strain energy rate EsDt     
					// Assumed constant dx and dy, since is located at this part of the code (i.e. end )
					EsDt = EsDt + 2*(shh_xx + shh_xy) * res_high*res_high*1; 

					// Exit if zeros are accessed, meaning that the calculation area is not properly filled
					//if (sii[m3]==0){fprintf(fp_log,"ERROR: Accessed a value in dissEs that was not filled !!!; m3=%d, gx(m1)=%e, gy(m2)=%e; while HR limits are m10=%d, m11=%d, m20=%d, m21=%d \n",m3,gx[m1],gy[m2],m10_hr,m11_hr,m20_hr,m21_hr); exit(0);}
				}
			}
		}
				
		// Final averages
		pra = pra / count;
		siia = siia / count;
		siia_nd = siia_nd / count;
		sxxa = sxxa / count;
		sxya = sxya / count;
		eiia = eiia / count;
		eiia_ne = eiia_ne / count;
		exxa_ne = exxa_ne / count;
		exya_ne = exya_ne / count;
		
		// Calculate effective area n m2 (since not all of defined had that rocktype)
		area = count*(res_high*res_high);
		
		// Integrate over time -> Strain energy Es
		Es1[na] = Es1[na] + EsDt*timestep;

        	if (count1 == 1)
        	{
                	sprintf(fileDissEs,"diss_es_%s_A%d.txt",exp_name,na);
                	FILE *fileexists = fopen(fileDissEs,"r");
                	if (fileexists)
                	{
                        fprintf(fp_log,">>> strain energy dissipation: loading file << %s >>  \n",fileDissEs); fflush(fp_log);
                	}
                else
                	{
                        fprintf(fp_log,">>> strain energy dissipation: new file << %s >> \n",fileDissEs); fflush(fp_log);
                        fld = fopen(fileDissEs,"w");
                        fprintf(fld,"Dimensions area (km) : %f %f %f %f \n",gx[m10_ca]/1e3,gx[m11_ca]/1e3,gy[m20_ca]/1e3,gy[m21_ca]/1e3);
			fprintf(fld," # t(Yr)          Es(N*m)        EsDt(N*m/s)       Pr(Pa)       siia(Pa)       siia_nd(Pa)       sxxa(Pa)      sxya(Pa)      eiia(s-1)      eiia_ne(s-1)      exxa_ne(s-1)      exya_ne(s-1)      area \n");
			fclose(fld);
                	}
        	}
		sprintf(fileDissEs,"diss_es_%s_A%d.txt",exp_name,na);
		fld = fopen(fileDissEs,"a");
		fprintf(fld,"%d  %f  %f  %f  %f  %f  %f  %f  %f  %e  %e  %e  %e  %f \n",count1,timesum05/3.15576e7,Es1[na],EsDt,pra,siia,siia_nd,sxxa,sxya,eiia,eiia_ne,exxa_ne,exya_ne,area);
		fclose(fld);
	} 
	/* ---- END: DISSIPATED STRAIN ENERGY RATE ---- */

	/* ---- START: THRUST INTERFACE ANALYSIS ---- */
	// Prepare output file 
	if (count1 == 1)
	{
		sprintf(fileThrust,"eachdt_thrust_%s_SF0.bin",exp_name);
		sprintf(fileThrustAbove,"eachdt_thrustabove_%s_SF0.bin",exp_name);
		sprintf(fileThrustBelow,"eachdt_thrustbelow_%s_SF0.bin",exp_name);
		// --- Get  GPS markers to track/follow ---
		// In postproc since need to select them once already properly subducting
		sprintf(fileTxtOutputGPS,"eachdt_gpsmarker_%s.txt",exp_name);
		sprintf(fileTxtMarkerOutput,"eachdt_markertrack_%s.txt",exp_name);
#if setup==12
		sprintf(fileTxtDPMarkerOutput,"eachdt_DP_markertrack_%s.txt",exp_name);
		sprintf(fileTxtUPMarkerOutput,"eachdt_UP_markertrack_%s.txt",exp_name);
#endif	
		// - Check if already have gps markers that need to be tracked -
		// If not at start of simulation than reload markers to be followed from previous txt file	
		FILE *fileexists = fopen(fileTxtOutputGPS,"r");
		if (fileexists)
		{
			fprintf(fp_log,"PAY ATTENTION : LOADING GPS MARKERS TO FOLLOW FROM %s = RIGHT FILE !?!? \n",fileTxtOutputGPS); fflush(fp_log);	
			
			// Note needs to be 'fl' for usage of ffscanf	
			fl = fopen(fileTxtOutputGPS,"rt");
			if (fl==NULL)
				{ fprintf(fp_log,"ERROR: Can not open gps file : %s \n ",fileTxtOutputGPS); fflush(fp_log); exit(0); }
			for (im=0;im<nr_gpsmarkers;im++)
				{ffscanf(); mgpstrack[im] = atoi(sa);}
			fclose(fl);
		}
		
		// - Select GPS markers to be followed -
		else
		{
			fprintf(fp_log,"SELECTING GPS MARKERS TO FOLLOW -> IS NO RELOAD OK ???? \n"); fflush(fp_log);
			
			//fprintf(fp_log," xgps[0] = %f, xgps[1] = %f, xgps[44] = %f \n",xgps[0],xgps[1],xgps[44]);
			for (im=0;im<nr_gpsmarkers;im++)
			{
				ygps[im] = 30e3;
				for (mm1=0;mm1<marknum;mm1++)
				{
					if (markt[mm1]>1 && markt[mm1]<8 && markx[mm1]>xgps[im] && markx[mm1]<(xgps[im]+1000) && marky[mm1]<ygps[im])
					{
						mgpstrack[im]=mm1;
						ygps[im] = marky[mm1];
					}	
				}
				//fprintf(fp_log," selected gps markers for im = %d : y = %f m, x = %f m, markt = %d, markernum = %d \n",im,marky[mgpstrack[im]],markx[mgpstrack[im]],markt[mgpstrack[im]],mgpstrack[im]); fflush(fp_log);
			}		
				
			// Write to be tracked marker numbers to file
			flgps = fopen(fileTxtOutputGPS,"wt");
			for (im=0;im<nr_gpsmarkers;im++)
				{fprintf(flgps,"%d ", mgpstrack[im]);}
			fprintf(flgps,"\n");
			fclose(flgps);
		}	
		
		// Check if all gps markers filled ; if not within three cells have markers, than code should fail !
		for (im=0;im<nr_gpsmarkers;im++)	  
		{
			if (mgpstrack[im]==0)
			{ 
				fprintf(fp_log,"  ERROR: Some GPS marker numbers (i.e. %d) to be tracked are not assigned !!!\n",im);
				fflush(fp_log);
				exit(1);
			}
		}
	}
	
	// --- Select thrust interface from maximum visco-plastic strainrate ---
	ybegin = (int)(15e3/res_high);
	yend   = (int)(120e3/res_high);
	xend   = (int)(m1_trench + 300e3/res_high); 	
	localdiprange = 8; // in nr grid cells
	
	check_bound(0,xnumx-1,xend,"postproc xend");
	check_bound(0,ynumy-localdiprange-5-1,yend,"postproc yend");
	
	// Determine sizes variables for binary writing
	szint=sizeof(int); 
	szdouble=sizeof(double);
	
	// Allocate memory thrust arrays and initialize to 0
	gx_th = (double*)calloc(xend-m1_trench,sizeof(double));
	gy_th = (double*)calloc(xend-m1_trench,sizeof(double));
	vpar_th = (double*)calloc(xend-m1_trench,sizeof(double));
	vslip_th = (double*)calloc(xend-m1_trench,sizeof(double));
	vslip_ne_th = (double*)calloc(xend-m1_trench,sizeof(double));
	vx_th = (double*)calloc(xend-m1_trench,sizeof(double));
	vy_th = (double*)calloc(xend-m1_trench,sizeof(double));
	maxstrain = (double*)calloc(xend-m1_trench,sizeof(double));
	localdip = (double*)calloc(xend-m1_trench,sizeof(double));
	sii_th = (double*)calloc(xend-m1_trench,sizeof(double));
	nu_th = (double*)calloc(xend-m1_trench,sizeof(double));
	ythrust = (int*)calloc(xend-m1_trench,sizeof(double));
	addythrust = (int*)calloc(xnumx,sizeof(int));
	
	// --- On thrust interface ---
	// Open binary file
	fldT = fopen(fileThrust,"ab"); 
	if (fldT==NULL)
		{ fprintf(fp_log,"Starting to write thrust data to new file! \n "); fflush(fp_log);}
	else
		{ fprintf(fp_log,"WARNING: Thrust data will now be appended!\n "); fflush(fp_log);}

        /* ============================================================== */
        /*              eulerian nodes on thrust interface                */
        /* ============================================================== */
	// Determine thrust interface in terms of nodes and coordinates 
	for(m1=m1_trench; m1<xend; m1++)
	{
		for(m2=ybegin; m2<yend; m2++)
		{
			m3=m1*ynumy+m2;
			m4=(m1+diff_m1)*yresol_vis+(m2+diff_m2);
			// Calculate visco-plastic strain rate; use normal, Non-Deviatoric stresses
			exx_ne = (sxx_nd[m3]/(2*nd[m3]));
			exy_ne = (sxy[m3]/(2*nu[m3]));
			eii_ne = sqrt(pow(exy_ne,2) + pow(exx_ne,2));
			if(ncomp[m4]==3 || ncomp[m4]==4 || ncomp[m4]==7 || ncomp[m4]==8 && eii_ne>maxstrain[m1-m1_trench])
			{
				if (m2>=ybegin && m2<=yend)
				{
					ythrust[m1-m1_trench] = m2;
					xythrust = m3;
					//fprintf(fp_log,"yhtrust[m1-m1trench] = %d \n", ythrust[m1-m1_trench]); getchar();
				}
				addythrust[m1] = m2;
				maxstrain[m1-m1_trench] = eii_ne;
				//printf("addythrust[m1]= %d \n",addythrust[m1]);
			}
		}
		// Form arrays that can be appended in binary file after this loop, so have them together for reading in 
		gx_th[m1-m1_trench] = gx[m1];
		gy_th[m1-m1_trench] = gy[ythrust[m1-m1_trench]];
		ybegin=(MAXV(ybegin,addythrust[m1])-1);
                //printf("ybegin = %d \n",ybegin);
	}

	/* ============================================================== */
        /*                  main thrust physical markers                  */
        /* ============================================================== */
	
	// - Check if already have physics markers that need to be tracked -
	if (count1==1)
	{
		// If not at start of simulation (t=0 s) than reload markers to be followed from previous txt file	
		FILE *fileexists = fopen(fileTxtMarkerOutput,"r");
		// - Select physics markers to be followed -
		// Load from previous file
		if (fileexists)
		{
			fprintf(fp_log,"PAY ATTENTION : LOADING PHYSICS MARKERS TO FOLLOW FROM %s = RIGHT FILE !?!? \n",fileTxtMarkerOutput); fflush(fp_log);	
			// Note needs to be 'fl' for usage with ffscanf	
			fl = fopen(fileTxtMarkerOutput,"rt");
			if (fl==NULL)
				{ fprintf(fp_log,"ERROR: Can not open physics marker file : %s \n ",fileTxtMarkerOutput); fflush(fp_log); exit(0); }
			for (im=0;im<nr_physmarkers;im++)
				{ ffscanf(); mphystrack[im] = atoi(sa); }
			fclose(fl);
		}
		// Select new
		else
		{
			fprintf(fp_log,"SELECTING PHYSICS MARKERS TO FOLLOW -> IS NO RELOAD OK ???? \n");

                        fflush(fp_log);
                        /* start from y: 15 km */
                        ylim_phys=15e3;
                        for (im=0;im<nr_physmarkers;im++)
                        {
                                xlim_phys=gx[m1_trench]+im*4.8e3;
                                m1 = m1serch(xlim_phys);
                                /*ylim_phys=(gy_th[m1-m1_trench]);*/
                                /*fprintf(fp_log,"ylim_phys starts from Y = %f \n",ylim_phys);*/
                                for (mm1=0;mm1<marknum;mm1++)
                                {
                                        /* constrain x & y limits */
                                        if (markx[mm1]>(xlim_phys-res_high) && markx[mm1]<(xlim_phys+res_high) &&  marky[mm1]>=(ylim_phys-res_high))
                                        {
                                                /* check rock type */
                                                if (markt[mm1]==3 || markt[mm1]==4 || markt[mm1]==7 || markt[mm1]==8)
                                                {
                                                mphystrack[im]=mm1;
                                                /*fprintf(fp_log,"mphystrack[im]= %f --- marky[mm1]= %f \n",mphystrack[im],marky[mm1]);*/
                                                break;
                                                }
                                        }
                                }
                                fprintf(fp_log," selected physics markers on main thrust for im = %d : y = %f m, x = %f m, markt = %d, markernum = %d \n",im,marky[mphystrack[im]],markx[mphystrack[im]],markt[mphystrack[im]],mphystrack[im]);
                        }
                        fflush(fp_log);

			/* Write to be tracked marker numbers to file */
			flphys = fopen(fileTxtMarkerOutput,"wt");
			for (im=0;im<nr_physmarkers;im++)
				{fprintf(flphys,"%d ", mphystrack[im]);}
			fprintf(flphys,"\n");
			fclose(flphys);
		}	
		
		// Check if all physics markers filled ; if not within three cells have markers, than code should fail !
		for (im=0;im<nr_physmarkers;im++)	  
		{
			if (mphystrack[im]==0)
			{ 
				fprintf(fp_log,"  ERROR: Some physics marker numbers (i.e. %d) to be tracked are not assigned !!!\n",im); fflush(fp_log);
				exit(1);
			}
		}
	}
	
	/* ============================================================== */
	/*	        downgoing plate physical markers 		  */
	/* ============================================================== */

#if setup==12	

	// - Check if already have physics markers that need to be tracked -
	if (count1==1)
	{
		// If not at start of simulation (t=0 s) than reload markers to be followed from previous txt file	
		FILE *fileexists = fopen(fileTxtDPMarkerOutput,"r");	
		// - Select physics markers to be followed -
		// Load from previous file
		if (fileexists)
		{
			fprintf(fp_log,"PAY ATTENTION : LOADING PHYSICS MARKERS +++downgoing plate+++ from %s = RIGHT FILE !?!? \n",fileTxtDPMarkerOutput); fflush(fp_log);
			//Note needs to be 'fl' for usage with ffscanf
			fl = fopen(fileTxtDPMarkerOutput,"rt");
	                if (fl==NULL)
                                { fprintf(fp_log,"ERROR: Can not open downgoing plate physics marker file : %s \n ",fileTxtDPMarkerOutput); fflush(fp_log); exit(0); }
                        for (im=0;im<nr_DP_physmarkers;im++)
                                { ffscanf(); DP_mphystrack[im] = atoi(sa); }
                        fclose(fl);
		}
		// Select new
		else
		{
			fprintf(fp_log,"SELECTING PHYSICS MARKERS TO FOLLOW -> IS NO RELOAD OK ???? \n");
			fflush(fp_log);
			
			for (im=0;im<nr_DP_physmarkers;im++)
			{
				xlim_DP_phys=gx[m1_trench]-im*4e3;	
				m1 = m1serch(xlim_DP_phys);
				ylim_DP_phys=35e3;
				for (mm1=0;mm1<marknum;mm1++)
				{
				if (markx[mm1]>(xlim_DP_phys-res_high) && markx[mm1]<(xlim_DP_phys+res_high) && marky[mm1]>ylim_DP_phys && marky[mm1]<(ylim_DP_phys+res_high))
					{
						DP_mphystrack[im]=mm1;
						break;
					}	
				}
				//fprintf(fp_log," selected physics markers in the downgoing plate for im = %d : y = %f m, x = %f m, markt = %d, markernum = %d \n",im,marky[DP_mphystrack[im]],markx[DP_mphystrack[im]],markt[DP_mphystrack[im]],DP_mphystrack[im]);
			}			
			fflush(fp_log);

			// Write to be tracked marker numbers to file
                        DP_flphys = fopen(fileTxtDPMarkerOutput,"wt");
                        for (im=0;im<nr_DP_physmarkers;im++)
                                {fprintf(DP_flphys,"%d ", DP_mphystrack[im]);}
                        fprintf(DP_flphys,"\n");
                        fclose(DP_flphys);
                }
		// Check if all physics markers filled ; if not within three cells have markers, than code should fail !
                for (im=0;im<nr_physmarkers;im++)
                {
                        if (mphystrack[im]==0)
                        {
                                fprintf(fp_log,"  ERROR: Some physics marker in the downgoing plate (i.e. %d) to be tracked are not assigned !!!\n",im); fflush(fp_log);
                                exit(1);
                        }
                }
        }

        /* - - - - - - - end*/


	/* ============================================================== */
	/*                  upper plate physical markers                  */
	/* ============================================================== */

	// - Check if already have physics markers that need to be tracked -
	if (count1==1)
        {
                // If not at start of simulation (t=0 s) than reload markers to be followed from previous txt file  	
                FILE *fileexists = fopen(fileTxtUPMarkerOutput,"r");
                // - Select physics markers to be followed -
                // Load from previous file
                if (fileexists)
                {
                        fprintf(fp_log,"PAY ATTENTION : LOADING PHYSICS MARKERS +++ upper plate +++ from %s = RIGHT FILE !?!? \n",fileTxtUPMarkerOutput); fflush(fp_log);
                        //Note needs to be 'fl' for usage with ffscanf
                        fl = fopen(fileTxtUPMarkerOutput,"rt");
                        if (fl==NULL)
                                { fprintf(fp_log,"ERROR: Can not open +++ upper plate +++ physics marker file : %s \n ",fileTxtUPMarkerOutput); fflush(fp_log); exit(0); }
                        for (im=0;im<nr_UP_physmarkers;im++)
                                { ffscanf(); UP_mphystrack[im] = atoi(sa); }
                        fclose(fl);
                }
                // Select new
                else
                {
                        fprintf(fp_log,"SELECTING PHYSICS MARKERS TO FOLLOW -> IS NO RELOAD OK ???? \n");
                        fflush(fp_log);

                        for (im=0;im<nr_UP_physmarkers;im++)
                        {
                                xlim_UP_phys=gx[m1_trench]+im*6e3;
                                m1 = m1serch(xlim_UP_phys);
                                ylim_UP_phys=35e3;
                                for (mm1=0;mm1<marknum;mm1++)
                                {
                                if (markx[mm1]>(xlim_UP_phys-res_high) && markx[mm1]<(xlim_UP_phys+res_high) && marky[mm1]>ylim_UP_phys && marky[mm1]<(ylim_UP_phys+res_high))
                                        {
                                                UP_mphystrack[im]=mm1;
                                                break;
                                        }
                                }
                                //fprintf(fp_log," selected physics markers in the +++ upper plate +++ numbers im = %d : y = %f m, x = %f m, markt = %d, markernum = %d \n",im,marky[UP_mphystrack[im]],markx[UP_mphystrack[im]],markt[UP_mphystrack[im]],UP_mphystrack[im]);
                        }
                        fflush(fp_log);
			
			// Write to be tracked marker numbers to file
			UP_flphys = fopen(fileTxtUPMarkerOutput,"wt");
			for (im=0;im<nr_UP_physmarkers;im++)
				{fprintf(UP_flphys,"%d ", UP_mphystrack[im]);}
			fprintf(UP_flphys,"\n");
			fclose(UP_flphys);
		}	
	
		// Check if all physics markers filled ; if not within three cells have markers, than code should fail !
		for (im=0;im<nr_physmarkers;im++)	  
		{
			if (DP_mphystrack[im]==0)
			{ 
				fprintf(fp_log,"  ERROR: Some physics marker numbers in the upper plate (i.e. %d) to be tracked are not assigned !!!\n",im); fflush(fp_log);
				exit(1);
			}
		}
	}

	/* - - - - - - - end*/

#endif

	// Gather nodal data on thrust interface
	for(m1=m1_trench+localdiprange; m1<xend-localdiprange; m1++)
	{
		m3=m1*ynumy+addythrust[m1];
		
		// Calculate fault parallel velocity using local dip
		if(addythrust[m1-localdiprange]!=0 && addythrust[m1+localdiprange]!=0)
		{
			check_bound(0,xend-m1_trench-1,m1-m1_trench+localdiprange,"postproc localdipcalc");
			check_bound(0,xend-m1_trench-1,m1-m1_trench-localdiprange,"postproc localdipcalc");
			localdip[m1-m1_trench] = atan((gy_th[m1-m1_trench+localdiprange]-gy_th[m1-m1_trench-localdiprange])/(gx_th[m1-m1_trench+localdiprange]-gx_th[m1-m1_trench-localdiprange])); 
			
			check_bound(0,xend-m1_trench-1,m1-m1_trench,"postproc vpar calc: arg localdip[]");
			vpar_th[m1-m1_trench] = vx[m3]*cos(localdip[m1-m1_trench]) + vy[m3]*sin(localdip[m1-m1_trench]);	
		}
		
		// Form arrays that can be appended in binary file after this loop, so have them together for reading in 
		vslip_th[m1-m1_trench] = 2*eii[m3]*res_high ;
		check_bound(0,xend-m1_trench-1,m1-m1_trench,"postproc maxstrain");
		vslip_ne_th[m1-m1_trench] = 2*maxstrain[m1-m1_trench]*res_high ;
		
		sii_th[m1-m1_trench] = sii[m3];
		nu_th[m1-m1_trench] = nu[m3];
		vx_th[m1-m1_trench] = vx[m3];
		vy_th[m1-m1_trench] = vy[m3];
	}
	
	// Append to binary files
	fwrite(&m1_trench,szint,1,fldT);
	fwrite(&xend,szint,1,fldT);
	fwrite(&timesum05,szdouble,1,fldT); 
	fwrite(gx_th,szdouble,xend-m1_trench,fldT);
	fwrite(gy_th,szdouble,xend-m1_trench,fldT);
	fwrite(vpar_th,szdouble,xend-m1_trench,fldT); 	// this one or vslip...	
	fwrite(vx_th,szdouble,xend-m1_trench,fldT); 	// tmp	
	fwrite(vy_th,szdouble,xend-m1_trench,fldT); 	// tmp 	
	fwrite(vslip_th,szdouble,xend-m1_trench,fldT);  // this one or vslip_ne... 
	fwrite(vslip_ne_th,szdouble,xend-m1_trench,fldT); // if not, than also maxstrain back to single value
	fwrite(sii_th,szdouble,xend-m1_trench,fldT); 		
	fwrite(nu_th,szdouble,xend-m1_trench,fldT); 		
	fclose(fldT);
	
	// free memory of arrays
	free(gy_th);
	free(vpar_th);
	free(vx_th);
	free(vy_th);
	free(vslip_th);
	free(vslip_ne_th);
	free(sii_th); 
	free(nu_th);
	free(maxstrain);  
	
	// --- 6.4 km above thrust interface ---
	// Equivalent to 1 cm in lab
	gy_ath = (double*)malloc( sizeof(double) * (xend-m1_trench) * 1 );
	vpar_ath = (double*)malloc( sizeof(double) * (xend-m1_trench) * 1 );
	sii_ath = (double*)malloc( sizeof(double) * (xend-m1_trench) * 1 );
	nu_ath = (double*)malloc( sizeof(double) * (xend-m1_trench) * 1 );
	
	// Initialise all entries in array to be 0 
	memset( gy_ath, 0, sizeof(double) * (xend-m1_trench) );
	memset( vpar_ath, 0, sizeof(double) * (xend-m1_trench) );
	memset( sii_ath, 0, sizeof(double) * (xend-m1_trench) );
	memset( nu_ath, 0, sizeof(double) * (xend-m1_trench) );
	
	// Open binary file
	fldT2 = fopen(fileThrustAbove,"ab");
	if (fldT2==NULL)
		{ fprintf(fp_log,"Starting to write thrust data to new file !\n"); fflush(fp_log);}
	else
		{ fprintf(fp_log,"WARNING: Thrust data will now be appended !\n"); fflush(fp_log);}
	
	// Gather data 1 cm above on thrust interface
	for(m1=m1_trench; m1<xend; m1++)
	{
		m2 = addythrust[m1]-above_thrust;
		m3=m1*ynumy+m2;
		gy_ath[m1-m1_trench] = gy[m2];
		// Form arrays that can be appended in binary file after this loop, so have them together for reading in 
		vpar_ath[m1-m1_trench] = vx[m3]*cos(localdip[m1-m1_trench]) + vy[m3]*sin(localdip[m1-m1_trench]);	
		sii_ath[m1-m1_trench] = sii[m3];
		nu_ath[m1-m1_trench] = nu[m3];
	}
	


	// Append to binary files
	fwrite(&m1_trench,szint,1,fldT2);
	fwrite(&xend,szint,1,fldT2);
	fwrite(&timesum05,szdouble,1,fldT2);
	fwrite(gx_th,szdouble,xend-m1_trench,fldT2);
	fwrite(gy_ath,szdouble,xend-m1_trench,fldT2);
	fwrite(vpar_ath,szdouble,xend-m1_trench,fldT2); 	// this one or vslip...	
	fwrite(sii_ath,szdouble,xend-m1_trench,fldT2); 		
	fwrite(nu_ath,szdouble,xend-m1_trench,fldT2); 		
	fclose(fldT2);

	free(gy_ath);
	free(vpar_ath);
	free(sii_ath); 
	free(nu_ath);
	
	// --- 6.4 km below thrust interface ---
	// Equivalent to 1 cm in lab
	gy_bth = (double*)malloc( sizeof(double) * (xend-m1_trench) * 1 );
	vpar_bth = (double*)malloc( sizeof(double) * (xend-m1_trench) * 1 );
	sii_bth = (double*)malloc( sizeof(double) * (xend-m1_trench) * 1 );
	nu_bth = (double*)malloc( sizeof(double) * (xend-m1_trench) * 1 );
	
	// Initialise all entries in array to be 0 
	memset( gy_bth, 0, sizeof(double) * (xend-m1_trench) );
	memset( vpar_bth, 0, sizeof(double) * (xend-m1_trench) );
	memset( sii_bth, 0, sizeof(double) * (xend-m1_trench) );
	memset( nu_bth, 0, sizeof(double) * (xend-m1_trench) );
	
	// Open binary file
	fldT3 = fopen(fileThrustBelow,"ab"); 
	if (fldT3==NULL)
		{ fprintf(fp_log,"Starting to write thrust data to new file ! \n "); fflush(fp_log);}
	else
		{ fprintf(fp_log,"WARNING: Thrust data will now be appended !\n "); fflush(fp_log);}
	// Gather data 1 cm above thrust interface
	for(m1=m1_trench; m1<xend; m1++)
	{
		m2 = addythrust[m1]+above_thrust;
		m3=m1*ynumy+m2;
		gy_bth[m1-m1_trench] = gy[m2];
		// Form arrays that can be appended in binary file after this loop, so have them together for reading in 
		vpar_bth[m1-m1_trench] = vx[m3]*cos(localdip[m1-m1_trench]) + vy[m3]*sin(localdip[m1-m1_trench]);	
		sii_bth[m1-m1_trench] = sii[m3];
		nu_bth[m1-m1_trench] = nu[m3];
	}	

	// Append to binary files
	fwrite(&m1_trench,szint,1,fldT3);
	fwrite(&xend,szint,1,fldT3);
	fwrite(&timesum05,szdouble,1,fldT3);
	fwrite(gx_th,szdouble,xend-m1_trench,fldT3);
	fwrite(gy_bth,szdouble,xend-m1_trench,fldT3);
	fwrite(vpar_bth,szdouble,xend-m1_trench,fldT3); 	// this one or vslip...	
	fwrite(sii_bth,szdouble,xend-m1_trench,fldT3); 		
	fwrite(nu_bth,szdouble,xend-m1_trench,fldT3); 		
	fclose(fldT3);
	
	free(gx_th);
	free(gy_bth);
	free(vpar_bth);
	free(sii_bth); 
	free(nu_bth);
	free(localdip);
	
	/* ---- END: THRUST INTERFACE ANALYSIS ---- */
	free(ncomp);
	
	/* --- GPS OUTPUT --- */
	flgps = fopen(fileTxtOutputGPS,"a");
	for (im=0; im<nr_gpsmarkers; im++)
		{ fprintf(flgps,"%e %.13e %.13e %e %e %e %e %e \n", timesum05, markx[mgpstrack[im]], marky[mgpstrack[im]], markp[mgpstrack[im]], markxx[mgpstrack[im]],markxy[mgpstrack[im]],markexx[mgpstrack[im]],markexy[mgpstrack[im]]); }
	fclose(flgps);
	
	/* --- Physics OUTPUT on main thrust --- */
	flphys = fopen(fileTxtMarkerOutput,"a");
	for (im=0; im<nr_physmarkers; im++)
		{ fprintf(flphys,"%e %.13e %.13e %e %e %e %e %e %e \n", timesum05, markx[mphystrack[im]], marky[mphystrack[im]], markp[mphystrack[im]], msbrit[mphystrack[im]], markxx[mphystrack[im]],markxy[mphystrack[im]],markexx[mphystrack[im]],markexy[mphystrack[im]]); }
	fclose(flphys);
#if setup==12
	/* --- Physics OUTPUT downgoing plate --- */
	DP_flphys = fopen(fileTxtDPMarkerOutput,"a");
	for (im=0; im<nr_DP_physmarkers; im++)
		{ fprintf(DP_flphys,"%e %.13e %.13e %e %e %e %e %e %e \n", timesum05, markx[DP_mphystrack[im]], marky[DP_mphystrack[im]], markp[DP_mphystrack[im]], msbrit[DP_mphystrack[im]], markxx[DP_mphystrack[im]],markxy[DP_mphystrack[im]],markexx[DP_mphystrack[im]],markexy[DP_mphystrack[im]]); }
	fclose(DP_flphys);
        /* --- Physics OUTPUT upper plate --- */
        UP_flphys = fopen(fileTxtUPMarkerOutput,"a");
        for (im=0; im<nr_UP_physmarkers; im++)
                { fprintf(UP_flphys,"%e %.13e %.13e %e %e %e %e %e %e \n", timesum05, markx[UP_mphystrack[im]], marky[UP_mphystrack[im]], markp[UP_mphystrack[im]], msbrit[UP_mphystrack[im]], markxx[UP_mphystrack[im]],markxy[UP_mphystrack[im]],markexx[UP_mphystrack[im]],markexy[UP_mphystrack[im]]); }
        fclose(UP_flphys);
#endif	
	fprintf(fp_log,"Print results each timestep to binary...OK! \n");
	fflush(fp_log);

	/* - - - check if Event Picking Algorithm file exists - - - */

	if (count1 == 1)
	{
		sprintf(filePicks,"pick_events_%s.txt",exp_name);
		FILE *fileexists = fopen(filePicks,"r");
		if (fileexists)
		{
			fprintf(fp_log,">>> pick events: loading file << %s >>  \n",filePicks); fflush(fp_log);
		}
		else
		{
			fprintf(fp_log,">>> pick events - new file: << %s >>  \n",filePicks); fflush(fp_log);
			flPK = fopen(filePicks,"w");
			fprintf(flPK,"rocktype    t(Yr)     X(m)     Y(m)    Vx(m/s)     Vy(m/s)    Pr(Pa)    sxx(Pa)     sxy(Pa)     exx(s-1)    exy(s-1)    eta(Pa s)    yield_stress(Pa)    slip_vel(m/s)  \n");
			fclose(flPK);
		}
	}

	/* - - - check if Plastic Yielding Counter file exists - - - */

	if (count1 == 1)
	{
		sprintf(fileCount,"plastic_yielding.txt");
		FILE *fileexists = fopen(fileCount,"r");
		if (fileexists)
		{
			fprintf(fp_log,">>> Plastic Yielding Counter: loading file << %s >>  \n",fileCount); fflush(fp_log);
		}
		else
		{
			fprintf(fp_log,">>> Plastic Yielding Counter: new file << %s >> \n",fileCount); fflush(fp_log);
			flCR = fopen(fileCount,"w");
			fprintf(flCR,">>>(1)Plastic-Yielding-Counter------(2)Seismic_Counter \n");
			fclose(flCR);
		}
	}

	/* ========================================================== */
	/* 		    Picking Event Algorithm		      */
	/* ========================================================== */

	/* VP_thresh = velocity threshold */
	
	/*       domain         */
	/*  X_11 -------- X_21  */
	/*   |		   |    */
	/*   |		   |    */
	/*  Y_12 -------- Y_22  */

	flPK = fopen(filePicks,"a");
	int plastic_counter = 0;
	int seismic_counter = 0;
	for (mm1=0;mm1<marknum;mm1++)
	{
		/* second invariant of stress Sii = sqrt(Sxx^2+Sxy^2);*/
		marksii[mm1] = sqrt(pow(markxx[mm1],2) + pow(markxy[mm1],2));
  	   
		if (markt[mm1]>1 && markt[mm1]!=7 && marke[mm1]>0 && markx[mm1]>X_11 && markx[mm1]<X_21 && marky[mm1]<Y_12)
		{
			/* plastic counter */	
			plastic_counter++;
			if (mvslip[mm1]>VP_thresh)
			{
				seismic_counter++;
				fprintf(flPK,"%d %.13e %.13e %.13e %e %e %e %e %e %e %e %e %e %e %e %e\n", markt[mm1], timesum05, markx[mm1], marky[mm1], markvx[mm1], markvy[mm1], markp[mm1], markxx[mm1], markxy[mm1], markexx[mm1], markexy[mm1], markv[mm1], msbrit[mm1], mvslip[mm1], msii_old[mm1],marksii[mm1]);
				fprintf(fp_log,"| mark-type=%d | X=%e | Y=%e | Vslip=%e | Syield=%e | Sii=%e \n",markt[mm1],markx[mm1],marky[mm1],mvslip[mm1],msbrit[mm1],marksii[mm1]); fflush(fp_log);
			}  
		}
	}
	fclose(flPK);


        flPK = fopen(filePicks,"a");
        for (mm1=0;mm1<marknum;mm1++)
        {
                /* second invariant of stress Sii = sqrt(Sxx^2+Sxy^2);*/
                marksii[mm1] = sqrt(pow(markxx[mm1],2) + pow(markxy[mm1],2));

                if (markt[mm1]==7 && marke[mm1]>0 && markx[mm1]>X_11 && markx[mm1]<X_21 && marky[mm1]<85e+3)
                {
                        /* plastic counter */
                        plastic_counter++;
                        if (mvslip[mm1]>5e-9)
                        {
                                seismic_counter++;
                                fprintf(flPK,"%d %.13e %.13e %.13e %e %e %e %e %e %e %e %e %e %e %e %e\n", markt[mm1], timesum05, markx[mm1], marky[mm1], markvx[mm1], markvy[mm1], markp[mm1], markxx[mm1], markxy[mm1], markexx[mm1], markexy[mm1], markv[mm1], msbrit[mm1], mvslip[mm1], msii_old[mm1],marksii[mm1]);
                                fprintf(fp_log,"| mark-type=%d | X=%e | Y=%e | Vslip=%e | Syield=%e | Sii=%e \n",markt[mm1],markx[mm1],marky[mm1],mvslip[mm1],msbrit[mm1],marksii[mm1]); fflush(fp_log);
                        }
                }
        }
        fclose(flPK);

	/* === Plastic Yielding Counter === */
	flCR = fopen(fileCount,"a");
	//fprintf(fp_log,">>> No Plastic Yielding Counter: %d",plastic_counter); fflush(fp_log);
	fprintf(flCR,"%d     %d   \n",plastic_counter,seismic_counter);
	fclose(flCR);
}
