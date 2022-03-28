/* MACROS --------------------------------------------- */
/* Min, Max, Abs */
#define MINV(a,b) ((a)<=(b)? (a):(b))
#define MAXV(a,b) ((a)>=(b)? (a):(b))
#define ABSV(a)   ((a)>=0? (a):-(a))
#define _TRUE_ 1
#define _FALSE_ 0
#if ((H5_VERS_MAJOR==1) && (H5_VERS_MINOR==8))
        #define H5_1_8
#elif( (H5_VERS_MAJOR_1==1 ) && (H5_VERS_MINOR==6) )
        #define H5_1_6
#endif

/* Main Sizes */
#define MAXMAT  	60000000		 /* General Matrix Size */
#define MAXVAL  	60000000 		 /* General Matrix Size */
#define MAXXLN  	1900 			 /* Max vertical line num in Grid */
#define MAXYLN  	401 			 /* Max horizontal line num in Grid */
#define MAXXMR   	8 			 /* Max marker in cell for X direction */
#define MAXYMR   	8 			 /* Max marker in cell for Y direction */
#define MAXTMR   	200			 /* Max markers types */
#define MAXPOS   	200			 /* Max Pos on Line num for wi[], wn[] buffers */
#define MAXFLN   	1000			 /* Max Output file names Num */
#define MAXNOD    	MAXXLN*MAXYLN  		 /* Max Total nodes num for all grid */
#define MAXCEL    	(MAXXLN-1)*(MAXYLN-1)  	 /* Max Total nodes num for all grid */
#define MAXPAR    	MAXXLN*MAXYLN*4   	 /* Max Total par num for all grid */
#define MAXMRK   	35000000 		 /* Max total marker Num */
#define MAXBON    	MAXPAR 		 	 /* Max Total bondary condition Equations for all grid */
#define setup 		12  			 /* Model setup options: LAB: 0=lab, ... ; NATURE: 10=nature dynw,JGR, 11=nature offev,GRL */
#define ls_t0 		1 			 /* For Large-Scale model: 1 = start from t=0 (geometry as defined in init.t3c, use mode_ls_to.t3c), 0 = start from prn defined at line in mode_ls_stm.t3c defined by file.t3c */

#if (setup<10)
	#define nr_gpsmarkers  	   11 /* nr gps markers in surface array */
	#define nr_physmarkers     14 /* nr physics markers in thrust array */
#elif (setup==10)
        #define nr_gpsmarkers      45 /* nr gps markers in surface array */
        #define nr_physmarkers     55 /* nr physics markers in thrust array */
#elif (setup==11)
        #define nr_gpsmarkers      45 /* nr gps markers in surface array */
        #define nr_physmarkers     55 /* nr physics markers in thrust array */
#elif (setup==12)
        #define nr_gpsmarkers      55 /* nr gps markers in surface array */
        #define nr_physmarkers     55 /* nr physics markers in thrust array */
        #define nr_DP_physmarkers  55 /* nr physics markers in downgoing plate array */
        #define nr_UP_physmarkers  55 /* nr physics markers in upper plate array */
#endif
/* End MACROS --------------------------------------------- */

/* ARRAYS ------------------------------------------------- */
/* Processing+Service Arrays for solutions */
/* val0[] - matrix contents */
/* fre0[] - free member for lines */
/* bufv[] - buffer for matrix organisation */
/* lin0[] - line numbers for matrix contents */
/* num0[] - pos numbers for line in  val0[], pos0[] */
/* pos0[] - first pos numbers for line in  val0[], pos0[] */
/* sol0[],sol1[] - solution buffers */
/* wn[],wi[] - Line Num,Koef buffer */
double sol0[MAXPAR],sol1[MAXPAR];
double tol0[MAXNOD],tol1[MAXNOD];
double val0[MAXMAT],fre0[MAXPAR],bufv[MAXPAR];
int lin0[MAXMAT],num0[MAXPAR],cur0[MAXPAR];
int pos0[MAXPAR];
int bufn[MAXPOS]; 
int wn[MAXPOS];
double wi[MAXPOS];

/* PARDISO SOLVER ARRAYS */
/* ia[] - first coloumn of each row */
/* ja[] - coloumn position of each koef in the respective row from left to right */
/* a[]  - value of each koef */
/* b[] - Right hand side */
/* x[] - solutions */
int ia[MAXPAR];             
int ja[MAXVAL];                  
double b[MAXPAR];           
double x[MAXPAR];           
double a[MAXVAL];                 

/* Nodes information */
/* gx[], gy[] - Coordinates of gridlines */
/* nd[] - Viscosity for SIGxx, SIGyy in Stokes equation, Pa*sek */
/* nu[] - Viscosity for SIGxy in Stokes equation, Pa*sek */
/* mu[] - Standard viskosity for node */
/* ep[] - Surface trace */
/* et[] - Adiabatic modulus */
/* gg[], gd[] - Shear modulus */
/* ro[] - Density, kg/m^3 */
/* tk[] - Temperature K */
/* cp[] - Heat capacity, J/kg */
/* kt[] - Thermal conductivity koef, Wt/m/K */
/* ht[] - Heat sources, Wt/m3 */
/* vx[],vy[] - Cur circle Vx,Vy,  m/sek */
/* tk0[] - Last circle TK */
/* wa0[], wa1[] - Water content */
/* td[] - Thermodynamic data base */
/* bondm[] - bondary position in bondn[],bondv[] for cur par (0) - no bondary */
/* bondn[] - PAR1+1,PAR2+1,PAR3+1 num in bond equat: CURPAR=CONST+KOEF1*PAR1 */
/* bondv[] - CONST,KOEF1,KOEF2,KOEF3 val in bond equat: CURPAR=CONST+KOEF1*PAR1 */
/* pr[] - pressure, Pa */
/* pr0[] - last cycle pressure, Pa */
/* pkf[] - koef for pressure aproximation */
/* exx[], exy[] - strain rates for cells (exx[]) and nodes (exy[]) */
/* esp[] - stores plastic strain rate (+=active, -=accumulated, non-active); Old Taras: rotation rate for nodes */
/* sxx[], sxy[] - stresses for cells (sxx[]) and nodes (sxy[]) */
/* sii[], sxx_nd[], sii_nd[] - second invariant of the (deviatoric) stress, nd = non-deviatoric */ 
/* eii[] - second invariant of the strain rate tensor */ 
/* sii[], sii_nd[] - second invariant of the (deviatoric) stress, nd = non-deviatoric */ 
/* sxxe[], sppe[], sxye[] - old elastic stresses, pressures for cells (sxxe[], sppe[]) and nodes (sxye[]) */
/* exxe[], exye[] - old strain rates from markers cells (exxe[] and nodes (exye[]) */
/* dro[] - Dro/Dt density changes with time */
/* drp[] - Dro/DP density changes */
/* mvx[], mvy[] - interpolated marker velocity in velocity points */
/* mrx[], mry[] - interpolated marker density in velocity points */
// sbrit(0)[] - Interpolated nodal yield strength (previous timestep)
// Start with capital is used for temporary storage during parallel run
double gx[MAXXLN],gy[MAXYLN];
double nu[MAXNOD],nd[MAXNOD],sxxe[MAXNOD],sppe[MAXNOD],sbritn[MAXNOD],sxye[MAXNOD],mu[MAXNOD],ep[MAXNOD],et[MAXNOD],gg[MAXNOD],gd[MAXNOD],ro[MAXNOD];
double nu0[MAXNOD],nd0[MAXNOD],sxxe0[MAXNOD],sppe0[MAXNOD],sbritn0[MAXNOD],sxye0[MAXNOD],exxe[MAXNOD],exxe0[MAXNOD],exye[MAXNOD],exye0[MAXNOD],ep0[MAXNOD],et0[MAXNOD],gg0[MAXNOD],gd0[MAXNOD],ro0[MAXNOD];
double tk[MAXNOD],cp[MAXNOD],kt[MAXNOD],ht[MAXNOD];
double cp0[MAXNOD],kt0[MAXNOD],ht0[MAXNOD];
double vx[MAXNOD],vy[MAXNOD];
double tk0[MAXNOD], tk1[MAXNOD], tk2[MAXNOD], tk3[MAXNOD];
double wa0[MAXNOD], wa1[MAXNOD];
double mvx[MAXNOD],mvy[MAXNOD],mrx[MAXNOD],mry[MAXNOD];
double mvx0[MAXNOD],mvy0[MAXNOD],mrx0[MAXNOD],mry0[MAXNOD];
double bondv[MAXBON][4];
double pr[MAXNOD],pkf[20];
long int bondm[MAXPAR],bondn[MAXBON][3];
double td[360][360][15][5];
double exx[MAXNOD],dro[MAXNOD],exy[MAXNOD],esp[MAXNOD];
double sxx[MAXNOD],drp[MAXNOD],sxy[MAXNOD];
double sxx_nd[MAXNOD],sii[MAXNOD],sii_nd[MAXNOD],eii[MAXNOD];
double sxx0[MAXNOD],sxy0[MAXNOD];
double dro0[MAXNOD],drp0[MAXNOD];

/* Markers and Rock information */
/* Markers information */
/* markx[], marky[] - X,Y of markers */
/* markk[] - temperature in K for marker */
/* markd[] - marker density, kg/m3 */
/* marke[] - finite strain for marker */
/* markp[] - pressure for marker */
// msbrit - yield stress Pa
/* markexx[], markexy[] - strain rate for marker */
/* markxx[], markxy[] - deviatoric stress components for marker */
/* markv[] - marker advected viscosity */
// markwa[] - pointer to indicater water presence Y(1)/N(0)
/* markvx[],markvy[] - marker velocity */
/* markt[] - rock type of markers */
/* Information for different rock types */
/* markim[] -  Immobility of marker type  Y(1)/No(0) */
/* marknu[], markdh[], markdv[], markss[] markmm[] -  Koef in ductile rheology Eq */
/* markll[], marka0[], marka1[], markb0[] markb1[], marke0[], marke1[], markf0[], markf1[] - Koef in brittle rheology Eq */
/* markgg[] - elastic shear modulus, Pa */
/* markn0[], markn1[], marks0[], marks1[] - viscosity and stress limits for individual rock types */
/* markro[], markaa[], markbb[], markcp[], markkt[], markkf[], markkv[], markht[] - ro,aro,bro,cp,kt, Tkt, Pkt, ht */
// msbrit - yield stress Pa
double markx[MAXMRK],marky[MAXMRK];
double markk[MAXMRK],markd[MAXMRK],marke[MAXMRK],markw[MAXMRK];
double markxx[MAXMRK],markv[MAXMRK],markxy[MAXMRK];
double markp[MAXMRK],markexx[MAXMRK],markexy[MAXMRK],msbrit[MAXMRK],msii_old[MAXMRK],marksii[MAXMRK];
double markvx[MAXMRK],markvy[MAXMRK],mvslip[MAXMRK];
char markt[MAXMRK];
double marknu[MAXTMR],markdh[MAXTMR],markdv[MAXTMR],markss[MAXTMR],markmm[MAXTMR];
double markgg[MAXTMR],markll[MAXTMR],marka0[MAXTMR],marka1[MAXTMR],markb0[MAXTMR],markb1[MAXTMR],marke0[MAXTMR],marke1[MAXTMR],markf0[MAXTMR],markf1[MAXTMR];
double markro[MAXTMR],markaa[MAXTMR],markbb[MAXTMR],markcp[MAXTMR],markkt[MAXTMR],markht[MAXTMR],markkf[MAXTMR],markkp[MAXTMR];
double markn0[MAXTMR],markn1[MAXTMR],marks0[MAXTMR],marks1[MAXTMR];
int markim[MAXTMR],markwa[MAXMRK];

/* Service, Connection Buffers */
/* vxy[] - Cur Vx,Vy */
/* eps[] - Cur Eps tenzors */
/* errbuf[] - Error buffer */
/* nunu[] - NUik buffer */
/* sa[] - input string buffer */
/* fl1in[] - input data file name */
/* fl1itp - input data file type */
/* fl1otp - output data file type */
/* fl0out[] - output data All file names */
/* fl0otp[] - output data files types */
/* fl0stp[],fl0cyc[] - output data file X,T,time steps cyc0max for Out files */
/* fl1out[] - output data Cur file name */
/* xn[], cn[] - high order FD cordinates and coefficients */
/* *fl - Stream for File Input/Output. Note: in ffscanf can only read fl, as hard-coded in that routine ! */
// n0 - krug; counter in prn generation cycle
// fnri - output nr for prn
double vxy[20],errbuf[20],nunu[9][5],eps[100];
char sa[2500];
char fl1in[50],fl1out[50];
int fl1itp,fl1otp;
char fl0out[MAXFLN][50];
int fl0otp[MAXFLN];
double fl0stp[MAXFLN][10];
int fl0cyc[MAXFLN];
double xn[1000],cn[1000][10];
FILE *fl,*fl1,*flfric;
int n0,fnri;

/* Arrays and variables used for hdf5 addition */
// mcomp pointers - Composition interpolated from marker on visualization grid, size dynamically determined in savehdf5.c */
// m..f pointers - followed marker arrays for storage in hdf5
// nm - number of markers that needs to be followed for picking algorithm (f is for fluid markers)
// follow[] - marker array to identify which markers to store in hdf5
int *mcomp,*mcf,*mcff;
double *markxf,*markyf,*markef,*markexxf,*markexyf,*markxxf,*markxyf,*markvf,*markpf,*markxff,*markyff,*markwff;
int start_cond,nm=0,nmf=0,nm_old,nmf_old;
char follow[MAXMRK];

/* End ARRAYS --------------------------------------------- */



/* VARIABLES -------------------------------------------------- */

/* Grid Parameters */
/* xnumx,ynumy - num of lines in grid in X,Y directions */
/* xnumx1,ynumy1 - num of cells in grid for X,Y directions */
/* mnumx,mnumy - num of markers in one cell for X,Y directions */
/* xsize,ysize - size of grid in X,Y directions, m */
/* GXKOEF,GYKOEF - Gravitation g for X,Y directions, m/sek^2 */
/* pinit - pressure on the upper boundary, Pa */
/* xstpx,ystpy - X,Y steps for grid */
/* kfx,kfy,kfxx,kfxy,kfyy,numdx,numdy - Koef for numeric differentiation */
/* cellnum - total Cell Num */
/* nodenum - total Node Num */
/* marknum - total Marker Num */
/* rocknum - rock types num */
/* bondnum - bondary condition equation num */
/* nodenum2 - nodenum*2 */
/* nodenum3 - nodenum*3 */
// n_glayer etc - nodal nr top of gelatine lower layer, nodal nr ok start gelatine, nodal nr end gelatine (backstop)
// mXS_hr - limits High Resolution area; X; 1=x, 2=y; S; 0=start, 1=end
// res_high - min. grid size, automatically from init.t3c
long int xnumx,ynumy,mnumx,mnumy;
long int xnumx1,ynumy1;
double xsize,ysize;
double GXKOEF,GYKOEF;
double pinit;
double xstpx,ystpy;
double kfx,kfy,kfxx,kfxy,kfyy,numdx,numdy,mardx,mardy;
long int cellnum,nodenum;
long int marknum,bondnum;
int rocknum=0;
long int nodenum2, nodenum3;
int n_glayer, nstart_gel, nend_gel;
int m10_hr,m11_hr,m20_hr,m21_hr;
double res_high;

/* Service parameters */
/* printmod - print information on the monitor Y(1)/N(0) */
/* intermod - order of interpolation from nodes to markers 0-linear >=1-Non linear */
/* intermod1 - order of interpolation from markers to nodes 0-linear >=1-Non linear */
/* densimod - mode of  density calculation: 0-constant, 1-PT-dependent */
/* outgrid - marker move out of grid Y(0)/N(1) Orthogonal Only (2) */
/* fl0num - number of otput file Names */
/* pos0cur - Cur pos Counter in val0[] lin0[] */
// pemod - mode to print marker picking algorithm information to hdf5-file every time step; 1=print, 0=do not print
long int pos0cur=0,printmod;
int fl0num,intermod,intermod1,outgrid=0,densimod,pemod;

/*  Bul ? Termodynamic database parameters <loadjxxx.c> */
/* savefluid - 1 = yes; save x,y,t of fluid markers in hdf5 after start_cond */
double tkmin,pbmin,tkstp,pbstp;
int tknum,pbnum;
double tkmin1,pbmin1,tkstp1,pbstp1;
int tknum1,pbnum1;
double tkpor,zmpor,vyfluid,vymelt,dmwamin,tdeep,zdeep,dxwater,dywater,deserp,dyserp;
double lambfld,teclmin,teclmax;
int savefluid;

/* Errosion/Sedimentation parameters */
/* erosmod - errosion/sedimentation Y(1)/N(0)  */
/* eroslev - errosion level, m */
/* eroscon - constant erosion rate, m/s */
/* eroskoe - increment of erosion rate with elevation, 1/s */
/* sedilev - sedimentation level, m */
/* sedicon - constant sedimentation rate, m/s */
/* sedikoe - increment of sedimentation rate with burial, 1/s */
/* sedimnum - Num of cycles of sedimentation */
/* sedimcyc - Num of cycles of sedimentation for One Layer */
/* waterlev - water/air boundary, m */
/* basalty - basalt melting depth, m */
/* dehydrmin, dehydrmax - serpentine dehydration min, max depth, m */
int erosmod;
double eroslev,eroscon,eroskoe;
double sedilev,sedicon,sedikoe;
double waterlev,basalty,dehydrmin,dehydrmax,slopemax;
int sedimnum=0,sedimcyc=3;

/* Motion parameters */
/* cyc0max - Num of circles for cur calculation */
/* maxxystep - max Distance change for one time step, m */
/* xelvismin - min viscoelasticity factor */
/* maxtkstep - max Themperature for one time step, K */
/* maxtmstep - max Time for one time step, sek */
/* timebond - time limit for internal boundary conditions from start */
/* timestep - time step, sek */
/* timestepe - time step for viscoelasticity, sek */
/* timesum - time from start, sek */
/* timedir - direction of motion in time:1 - forward, -1 - backward */
/* movemod - solve Stoks+Continuity equations Y(1)/N(0) */
/* tempmod - solve Heat Transfer equation Y(1)/N(0) */
/* markmod - move markers Y(1)/N(0) */
/* ratemod - reset velosity and pressure Y(1)/N(0) */
/* gridmod - recalc grid parameters Y(1)/N(0) */
/* inertyn - inertia in the equations Y(1)/N(0) */
int cyc0max,movemod,tempmod,markmod,ratemod,gridmod;
double maxxystep,xelvismin,maxtkstep,maxtmstep,timestep,timestepe,timesum,timedir;
double timebond,inertyn;
/* General iteration parameters in <gaus.c> */
/* ckoef - Cur Koef for Zeydel method of iteration */
double ckoef;

/* V parameters in <move.c> */
/* DIVVMIN,STOKSMIN - Min Valid absolut Err val for Contin,Stokes Eq */
/* stoksmod - dNu/dT, dNu/dP Y(0)/N(1) in Stokes calc */
/* viscmod - Viscoelasticity Y(1,2)/N(0) */
/* stoksfd - Order of FD for EPSxx,EPSyy,EPSxy calculation */
/* nukoef - No function */
/* nubeg,nuend,nucontr - Min, Max, Max/Min limits of viscozity */
/* hidry - max depth with hidrostatic pressure of pore fluid */
/* hidrl - brittle weackening factor for hidrostatic pressure of pore fluid */
/* strmin, strmax - min max values of stress Pa allowed */
/* stredif - numerical stress relaxation coefficient */
double DIVVMIN,STOKSMIN;
int stoksmod,stoksfd;
double viscmod;
double nubeg,nuend,nucontr,hidry,hidrl,strmin,strmax,stredif;

/* T parameters in <heat.c> */
/* HEATMIN - Min Valid absolut Err val for Heat Equation */
/* heatmod - dK/dX account Y(1)/N(0) */
/* heatfd - Order of FD for Qx, Qy calculation */
/* heatdif - Numerical diffusion coefficient */
/* frictyn - Viscouse friction heat Y(1)/N(0) */
/* adiabyn - adiabatic heat calculation: N(0)/Y(1) */
double HEATMIN;
double heatdif;
int heatmod,heatfd,frictyn,adiabyn;
int collision;
double timecc0,timecc1,vel_cc1,vel_cc2,vel_cc,vel_cc0;
// Global variables added for STM in the lab
// count1 	- counter for number of output text file lines, i.e. number of time steps since start outputting
// restart	- switch to determine if restart and need to load prior markers followed 
// debugmod - 0 do not print output useful for debugging, 1 do print it
// file...	- filenames used for output data each timestep
// vb_layer,above_fric - Y nodal point locators for each timestep output
// metis_reorder- switch to determine if can do metis reordering of global matrix; not if just after restart as potentially unstable
// veldepfric - switch for turning rate-dependent friction on and off
// gelx0 - x-coordinate of point 0 of the gel-wedge as defined in the init.t3c file
// gely0 - y-coordinate of point 0 of the gel-wedge as defined in the init.t3c file
// w_height - height of the wedge
// d_bstop - distance from the backstop for marker m7 (track markers:figure 7 in vanDinther (2013))
// startgps - relative distance gps marker selection start w.r.t. gelx0
// dxgps - distance between gps markers
// mus - Static friction coefficent
// mgamma - amount of frictional weakening allowed (gamma)
// mvc - Characteristic Velocity, determining at which velocity friction is reduced by a half
// .._: mt - MegaThrust; vw = Velocity-Weakening, vs = Velocity-Strengthening
// start_sez->end_downdip- values determining location seismogenic zone
// slab_dip(_deg) - slab dip for picking which (surface) markers to follow
// vpush	- push velocity for calculating actual slip, assuming motion subducting plate analogue is constant
// vtresh	- coseismic velocity threshold for event selection:
// event	- switch to determine if during event (and event_nr) or not for output accumulating slip -acc_slip's-
// mphystrack	- array with (near) interface marker tracked
// mgpstrack	- array with surface marker tracked
// sbrit_ave 	- strength averaged over markers within seismogenic zone
// count_sezm 	- counter to calculate nr of markers within seismogenic zone
// fp_log	- where to print output to directly (useful during parallel calculations)
int count1=0, restart=0, debugmod;
char fileTxtOutput[50], fileTxtMarkerOutput[50],fileTxtDPMarkerOutput[50], fileTxtUPMarkerOutput[50], fileTxtOutputS[50], fileTxtOutputE[50], fileTxtOutputP[50],fileTxtOutputQ[50], fileTxtOutputSxx[50], fileTxtOutputSxy[50], fileTxtOutputDRg[50],fileTxtOutputDRfd[50], fileTxtOutputDRfs[50], fileTxtOutputGPS[50], file2open[50];
int vb_layer,above_fric,metis_reorder=0,veldepfric;
double gelx0, gely0, w_height, d_bstop, startgps, dxgps, slab_dip_deg;
double mus_mtvw,mus_vs,mgamma_vw,mgamma_vs,mvc_vw,mvc_vs,mus_vs;
double start_sez, end_sez, before_trench, half_range, start_updip, end_updip, start_downdip, end_downdip, slab_dip, vpush, vtresh, event, *acc_slip_gel,*acc_slip_fbld, *acc_slip_fbls;
int mphystrack[nr_physmarkers], mgpstrack[nr_gpsmarkers], event_nr;
#if (setup==12)
int DP_mphystrack[nr_DP_physmarkers], UP_mphystrack[nr_UP_physmarkers];
#endif
double sbrit_ave,count_sezm;
const char log_filename[] = "STM.out";
FILE *fp_log; 

// Additional STM global variables for large-scale
// meltmod - switch melting on(1) or off(0)
// above_thrust - number of nodes above thrust where velocity is analyzed
// tk_updipsez0/1 - up- and downdip limits of seismogenic zone in terms of temperature
// dX/Ygps - arrays with x and y coordinates of markers tracking gps/surface displacements
// shift_km - value in km by which grid should be shifted to conform high resolution area to trench
// mus_ini - static friction coefficient during INItial timesteps (before STM)
// lambfld_ini - porefluid pressure ratio during INItial timesteps (before STM)
// lambfld_stm - porefluid pressure ratio during STM
// file... - filenames for frequent output
int above_thrust, meltmod;
double tk_updipsez0,tk_updipsez1,xgps[nr_gpsmarkers],ygps[nr_gpsmarkers],shift_km,mus_ini,lambfld_ini,lambfld_stm;
char fileDissEsBasic[50],fileThrust[50],fileThrustAbove[50],fileThrustBelow[50],fileTxtOutputFric[50],filePicks[50],fileCount[50],exp_name[20];
double X_11,X_21,Y_12,VP_thresh;

/* End VARIABLES -------------------------------------------------- */



/* FUNCTONS PROTOTYPES ------------------------------------------- */
// omp at end of a function name means that this routine has been adapted for Open MP parallelization

/* <load.c> FILE */
/* ffscanf() - Load single word from file fl without empty lines */
/* ffscanf1() - Load single word from file fl1 without empty lines. Note these 2 routines are exactly the same, except for the filename they read from. */
/* loadconf() - Load configuration from mode.t3c */
/* loader() - Load information from data file */
/* saver() - Save information to data file */
/* gridcheck() - Calc,Check parameters of Grid */
/* checkbound() - check for boundaries array */
void ffscanf();
void ffscanf1();
int loadconf();
void loader();
void saver(int, int);
void gridcheck();
void check_bound(int ,int , int , const char *);

/* <gaus.c> FILE */
/* gausmat3() - Solve system of linear equation by economic frontal Gauss method */
/* gausmat4() - Solve system of linear equation by Pardiso solver */
int gausmat3(int, long int, long int);
int gausmat4(int, long int, long int);

/* <move.c> FILE */
/* viterate() - General subroutine for Vx,Vy,P Calc by Iterativ method */
/* xstokserr() -  Right part or  Err in X Stokes Equat for cur node calc */
/* ystokserr() -  Right part or  Err in Y Stokes Equat for cur node calc */
/* conterr() -  Right part or  Err in Contin Equat for cur node calc */
/* xbonderr() -  Right part or  Err in Boundary vX Equat for cur node calc */
/* ybonderr() -  Right part or  Err in Boundary vY Equat for cur node calc */
/* pbonderr() -  Right part or  Err in Boundary P Equat for cur cell calc */
/* sxxcalc() etc. - Value or add EPS and SIG equations */
/* fdweight() - weight for FD calculation after Fornberg (1996) */
void viterateomp(int);
double xstokserr(long int, long int, int);
double ystokserr(long int, long int, int);
double conterr(long int, long int, int);
double xbonderr(long int, int);
double ybonderr(long int, int);
double pbonderr(long int, int);
double sxxcalc(long int, long int, double);
double sxycalc(long int, long int, double);
void fdweight(int, int, double);

/* <heat.c> FILE */
/* titerate() -  Temperature recalc after time step */
/* heatserr() -  Right+Right part or  Err in Heat Equat for cur node calc */
/* tbonderr() -  Right part or  Err in Boundary T Equat for cur node calc */
/* tkrecalc() -  Calc average temperature after new marker position */
/* qxcalc(),qycalc() - Coefficients or value for Qx,Qy  Equations */ 
void titerateomp(int);
double heaterr(long int, long int, int);
double tbonderr(long int, int);
void tkrecalc();
double qxcalc(long int, long int, double);
double qycalc(long int, long int, double);

/* <mark.c> FILE */
/* movemark() - move markers by Runge-Kutta method */
/* ronurecalc() - recalc ro[],nu[] etc after new marker position */
/* dencalc() -  Calc density for given P,T after ro equation */
/* tdbasecalc() - computing TD variables from the database */
/* antigor() - Mantle transformation to Antigorite */
/* viscalc() -  Calc viscosity for given P,T after rheol equation */
/* hydration2() - Hydration front progress recalc */
/* melting(), meltpart(), meltpart1() - melting of rocks account */
/* erosion() - Erosion Surface Tracing */
/* erosmark() - Erosion/Sedimentation Function for markers */
/* m1serch(), m2serch() - serch of nearest upper-left node of cell for  current location */
/* allinteri() - Vx,Vy, EPS calc for marker by interpolation */
/* allintert() - T calc for marker by interpolation */
/* allinterp() - P calc for marker by interpolation */
/* allinterd() - Sxx,Syy,Sxy calc for markers */
/* allinters() - Vx,Vy,EPS*SIG calc for marker by interpolation */
/* nodewt() - Weights for horisontal and vertical nodes calculation for marker interpolation */ 
void movemarkomp();
void ronurecalcomp();
void dencalcomp(double, double, double, double, int, double *, double *, double *);
void antigoromp(double, double, double, double, long int, long int, char *);
void viscalcomp(double, double, double, double, double, double, double, double, double, double, double, double *, double *, double *, long int, int, int, long int,double *, double *);
double hydration2omp();
void tdbasecalcomp(double, double, double, double, int, long int, long int, double *, double *, double *, double *, double *, double *, double *, double *);
void erosmarkomp(long int, int, long int, double, double, char *, double *, double *);
void meltingomp(double, double, long int, int, char *, double *, double *, double *);
void meltpart1omp(double, double, int, double *, double *);
void meltpartomp(double, double, double, double, long int, int, double *,double *, double *, double *, double *, double *, double *, double *, double *);
void erosion();
long int m1serch(double);
long int m2serch(double);
void nodewt(long int, long int, long int, long int, double, double, int, int);
void allinteriomp(double, double, long int, long int, double *, double *, double *, double *, double *);
void allinterdomp(double, double,long int, long int, double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *, double *, double *, double *);
void allintertomp(double, double, long int, long int, double *,double *);
double allinterpomp(double, double, long int, long int);
void allinters(double, double);
void allinteri(double, double);

/* <hdf5.c> FILE */
// Routines to generate hdf5 files
void create_output(const char*);
void AddGroup_to_hdf5(const char*, const char*);
void AddFieldToGroup_generic(int, const char*, const char*, const char*, char, int, void*, int);
int create_hdf5(int,char*,int);
int create_hdf5_markerprop(int,char*,int);
int interpoler(int*, int, int, int);
double spncalc(long int, long int);


/* <spec*.c> FILE */
// STM add-ons
// flyinput - input parameter that can overwrite on the fly, after loading prn data
// set_stmsw - sets global variables for large-scale stm depending on the time step
// postproc - save post-processing data each time step to text file
// preppostanalysis - check wether natural model far enough to start analysis
void flyinput_lab(); void load_flyinput_lab(); void flyinput_dynw(); 
void set_stmsw();
void postproc_lab(); void postproc_dynw(int *,double *,int);
void preppostanalysis(double **,int *);

/* End FUNCTONS PROTOTYPES ------------------------------------------- */


