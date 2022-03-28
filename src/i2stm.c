#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <mkl.h>
#include <hdf5.h>
#include <omp.h>
#include <ctype.h>
#include <time.h>

#include"headstm.c"
#include"loadstm.c"
#include"movestm.c"
#include"markstm.c"
#include"gausstm.c"
#include"heatstm.c"
#include"hdf5stm.c"
#include"speclsstm.c"

//=======  lab code =========
//#include"speclabstm.c"
//#include"extrainputlabstm.c"

/* Solve differential equations by MARKER-IN-CELL + FINITE-DIFFERENCES method */
int main()
{
	// Timers for profiling
	clock_t start=0, cur_t=0, t_init=0, t_load=0, t_vpit=0, t_maxv=0, t_tite=0, t_move=0, t_ronu=0, t_save=0, t_post=0;

	// Decleration
	int n1,f0,n3,fln3,pos,pos1,nr_areas;
	long int pos0cur0,m1,m2,m3;
	char keys[] = "1234567890";
	char numbers[] = ".";
	char txtprn[50],txthdf5[50],txthdf5_marker[50],txthdf5_fluid[50];
	double *Es1;
	
	// Open log file for all output, so can flush and always have output written
	fp_log = stdout;

	if (printmod){fprintf(fp_log,"Starting... \n"); fflush(fp_log);}

	/* Load configuration from mode.t3c */
	fln3=loadconf()+1;
	if (printmod){fprintf(fp_log,"Loaded configuration mode.t3c. \n"); fflush(fp_log);}

	// Make filename for additional output: cut experiment name out of output file name
	for (n1=0;n1<50;n1++) fl1out[n1]=fl0out[1][n1];
	pos = strcspn (fl1out,keys);
	strncpy (exp_name,fl1out,pos);
	exp_name[pos]='\0';
	puts (exp_name);

	// Form new file name for output every timestep
#if setup<10
	sprintf(fileTxtOutput,"eachdt_%s.txt",exp_name);
	sprintf(fileTxtMarkerOutput,"eachdt_markertrack_%s.txt",exp_name);
	sprintf(fileTxtOutputGPS,"eachdt_gpsmarker_%s.txt",exp_name);
	sprintf(fileTxtOutputS,"eachdt_sii_%s.txt",exp_name);
	sprintf(fileTxtOutputQ,"eachdt_qy_%s.txt",exp_name);
	sprintf(fileTxtOutputP,"eachdt_pr_%s.txt",exp_name);
	sprintf(fileTxtOutputE,"eachdt_eii_%s.txt",exp_name);
#endif
		
	/* Load data from prn file */
	cur_t = (double)omp_get_wtime();
	loader();
	t_load = (double)(omp_get_wtime() - cur_t);
	
	// Initialize counters
	event = 0; event_nr = 0;


	/* --- Output file cycle --- */
	for (f0=fln3;f0<fl0num;f0++)
	{
		/* Reload current output file name, type */
		for (n1=0;n1<50;n1++) fl1out[n1]=fl0out[f0][n1];
		fl1otp=fl0otp[f0];

		// Reload cyc0max_maxxystep_maxtkstep_maxtmstep_maxdsstep for this output cycle 
		cyc0max=fl0cyc[f0];
		maxxystep=fl0stp[f0][0];
		maxtkstep=fl0stp[f0][1];
		maxtmstep=fl0stp[f0][2];
		xelvismin=fl0stp[f0][3];

		// Set switches for inertia, post-processing, and strength depending on time step 
		// Allows to run in once from one mode.t3c up to point of analysis start
#if setup > 9
		set_stmsw();
#endif	
		
		// If trouble with PARDISO stability: use reordering of METIS=3 after 10 timesteps from (re)start to speed up again
		// METIS=3 is not stable on monch
		//if (f0>=fln3+10){metis_reorder = 1;}
		//else {metis_reorder = 0;}

		/* --- Timestepping loop within each output file --- */
		for (n0=0;n0<cyc0max;n0++)
		{
			if (printmod) fprintf(fp_log,"\n! Start on FILE %s  KRUG %d ! \n",fl1out,n0+1); fflush(fp_log);
			
			// Set average sez counter to zero 
			sbrit_ave = count_sezm = 0;

			// Dynamic Rupture-like output this model: with nr so if restart print to new file
			if (setup<10 && count1==0)
			{
				sprintf(fileTxtOutputDRg,"dr_gel_%s%i.txt",exp_name,f0);
				sprintf(fileTxtOutputDRfs,"dr_fbls_%s%i.txt",exp_name,f0);
				sprintf(fileTxtOutputDRfd,"dr_fbld_%s%i.txt",exp_name,f0);
			}
			else if (setup>9 && n0==1)
			{
				sprintf(fileTxtOutputFric,"fricprop_channel_%s_fr%d.txt",exp_name,f0+1);
			}

			/* *** Step 1: Set initial time step *** */
			timestep=maxtmstep;			// Displacement timestep
			timestepe=xelvismin;		// Computational elastic timestep
			t_init=(double)omp_get_wtime();
			if (printmod) fprintf(fp_log,"\n !!! MAX VALID TIME STEP FOR CYCLE %e YEARS !!! \n",timestep/3.15576e+7); fflush(fp_log);

			/* *** Step 2-5: vX,vY recalc after Stokes+Continuity equation *** */
			if(movemod)
			{
				if (printmod) fprintf(fp_log,"\n EPS, SIG, P, VX, VY CALCULATION...\n"); fflush(fp_log);
				cur_t = (double)omp_get_wtime();
				viterateomp(n0);
				t_vpit = (double)(omp_get_wtime() - cur_t);
				if (printmod) fprintf(fp_log,"EPS, SIG, P, VX, VY  OK!\n"); fflush(fp_log);
			}

			/* *** Step 6-9: Calculate temperature following heat transport equation *** */
			// Excluded in lab setup due to tempmod in mode.t3c
			if(timedir>0 && tempmod && timestep)
			{
				if (printmod) fprintf(fp_log,"\n TEMPERATURE CALCULATION...\n"); fflush(fp_log);
				cur_t = (double)omp_get_wtime();
				titerateomp(n0);
				t_tite = (double)(omp_get_wtime() - cur_t);
				if (printmod) fprintf(fp_log,"TEMPERATURE OK!\n"); fflush(fp_log);
			}

			/* *** Step 10: Move marker *** */
			if(markmod && timestep)
			{
				/* Time step for markers definition */
				if(timedir<0 && timestep>0) timestep=-timestep;
				if (printmod) fprintf(fp_log,"\n CURRENT VALID TIME STEP %e YEARS IN CYCLE\n",timestep/3.15576e+7); 
				if (printmod) fprintf(fp_log,"MOVE MARKERS..."); fflush(fp_log);
				cur_t = (double)omp_get_wtime();
				movemarkomp();
				t_move = (double)(omp_get_wtime() - cur_t);
				if (printmod) fprintf(fp_log,"MARKERS MOVED OK!\n"); fflush(fp_log);
			}

			/* *** Finish up *** */

			/* Reset Vx, Vy, P values */
			if(ratemod)
			{
				if (printmod) fprintf(fp_log,"VX, VY, P RESET ..."); fflush(fp_log);
				
				/* Vx, Vy Reset Cycle */
				for (m1=0;m1<nodenum;m1++)
					{vx[m1]=vy[m1]=0;}

				/* Pressure in cells Reset Cycle */
				for (m1=0;m1<xnumx1;m1++)
					for (m2=0;m2<ynumy1;m2++)
				{
					// Pos in pr[]
					// ATAT TARAS Why use ynumy1 instead of ynumy? ynumy1=ynumy-1, so one less should go wrong, no?
					m3=m1*ynumy1+m2;
					// Recalc P
					pr[m3]=pinit+((double)(m2)+0.5)*ystpy*GYKOEF*pkf[0];
				}
				if (printmod) fprintf(fp_log,"VX, VY, P OK!"); fflush(fp_log);
			}
		
			/* Post-processing: save some information every timestep */
			cur_t = (double)omp_get_wtime();
#if setup < 10
			start_cond = 1;
			postproc_lab();
#elif setup > 9
			// Start condition defined by size of timestep in i2.c
			if (start_cond==1)
			{
				if (count1==0)
					{ preppostanalysis(&Es1,&nr_areas); }
				postproc_dynw(&nr_areas,Es1,f0);
			}
#endif 
			t_post = (double)(omp_get_wtime() - cur_t);

			/* Save marker results to hdf5 every timestep for picking algorithm */
			if(pemod==1 && start_cond==1 && fl1otp==2)
			{
				sprintf(txthdf5_marker,"%s%1.3i_dt%1.2i_marker",exp_name,f0,n0);
				create_hdf5_markerprop(_TRUE_, txthdf5_marker, n0);
			}

			/* Increase time (1 year = 3.15576*10^7 s) */
			timesum+=timestep;
			if (printmod) fprintf(fp_log,"\n %e YEARS IN CYCLE     %e YEARS FROM START\n\n",timestep/3.15576e+7,timesum/3.15576e+7); fflush(fp_log);

			/* Print times for profiling */
			if (printmod) 
			{
				fprintf(fp_log,"--------------------\n");
				fprintf(fp_log,"Wall time:       %3.3f sec\n", (double)(omp_get_wtime() - t_init));
				fprintf(fp_log,"--------------------\n");
				//fprintf(fp_log,"loader took:     %3.3f sec\n", (double)t_load);
				//fprintf(fp_log,"saver took:      %3.3f sec\n", (double)t_save);
				//fprintf(fp_log,"--------------------\n");
				//fprintf(fp_log,"ronurecalc took: %3.3f sec\n", (double)t_ronu);
				fprintf(fp_log,"movemark took:   %3.3f sec\n", (double)t_move);
				fprintf(fp_log,"--------------------\n");
				fprintf(fp_log,"vpiterate took:  %3.3f sec\n", (double)t_vpit);
				//fprintf(fp_log,"titerate took:   %3.3f sec\n", (double)t_tite);
				fprintf(fp_log,"--------------------\n");
				fprintf(fp_log,"postproc took:   %3.3f sec\n", (double)t_post);
				fprintf(fp_log,"--------------------\n");
				fflush(fp_log);
			}
		}
		/* End time stepping loop */

		/* --- Optimize output --- */
		cur_t = (double)omp_get_wtime();
		if(fl1otp==2)
		{
			if (printmod) fprintf(fp_log,"Optimize output with producing hdf5 and prn removing\n"); fflush(fp_log);
			
			// Settings for removing prn files to save storage space
			int keepeveryxprn = 25;
			int prnoffset = 3;							// keep last 'prnoffset' prn's for possible restarts, remove those before
			
			// Make filenames
			int f0rm = f0-prnoffset;
			sprintf(txtprn,"%s%i.prn",exp_name,f0rm+1);
			sprintf(txthdf5,"%s%1.3i",exp_name,f0+1);
			sprintf(txthdf5_fluid,"%s%1.3i_fluid",exp_name,f0+1);
	
			// Remove previous .prn file from the directory, unless each 25th so can later restart from these (or last 5)
			if (((float)f0rm/keepeveryxprn)==roundf(f0rm/keepeveryxprn))
				{fprintf(fp_log,"Keep prn-file %s for possible later restarts. \n",txtprn); fflush(fp_log);}
			else
			{
				FILE *fileexists = fopen(txtprn,"r"); 
				if (fileexists) 
				{
					fclose(fileexists);
					int an = remove(txtprn);
					if( an != 0 )
						{fprintf(fp_log,"WARNING: prn-file %s can not be deleted\n",txtprn); fflush(fp_log);}
					else
						{fprintf(fp_log,"Prn-file %s successfully deleted\n",txtprn); fflush(fp_log);}
				}
				else if (fileexists == NULL) 
					{fprintf(fp_log,"WARNING: Prn-file %s to be removed can not be opened! \n",txtprn); fflush(fp_log);}
			}

			// Print results to hdf5
			create_hdf5(_TRUE_, txthdf5, n0);
			if (start_cond==1 && savefluid==1) {create_hdf5_fluid(_TRUE_, txthdf5_fluid, n0);}
		}

        // Print results to prn (contains all information for restarting)
		saver(f0+1,n0-1);
		t_save = (double)(omp_get_wtime() - cur_t);

		if (printmod) fprintf(fp_log,"End of 1 frame round as defined in mode.t3c.\n"); fflush(fp_log);
		/* End Output file Names Cycle */
	}

	/* End Program */
	if (printmod) fprintf(fp_log,"Your STM simulation ended successfully!\n"); fflush(fp_log);
	return 0;
}
/* End solve differential equations by MARKER-IN-CELL + FINITE-DIFFERENCES method */

