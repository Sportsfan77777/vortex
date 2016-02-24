#include "mp.h"

boolean         TimeToWrite, Restart = NO, OpenInner = NO;
int             TimeStep = 0, begin_i = 0, NbRestart = 0, verbose = NO;
static int      InnerOutputCounter=0, StillWriteOneOutput;

extern boolean  Corotating;
real            ScalingFactor = 1.0;
real            LostMass, LostLabel;


int
main(argc, argv)
int argc;
char *argv[];
{
  PolarGrid   *gas_density, *gas_v_rad, *gas_v_theta, *gas_label;
  PolarGrid1D   *gas_density_1D, *gas_v_rad_1D, *gas_v_theta_1D, *gas_label_1D;
  int          i;
  boolean      disable = NO, TimeInfo = NO;
  boolean      TestRadii = NO, Profiling = NO;
  TimeProcess  t_Hydro;
  char         ParameterFile[256];
  PlanetarySystem *sys;

  MPI_Init (&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &CPU_Rank);
  MPI_Comm_size (MPI_COMM_WORLD, &CPU_Number);
  CPU_Master = (CPU_Rank == 0 ? 1 : 0);
  IAmTheFirst = (CPU_Rank == 0 ? 1 : 0);
  IAmTheLast  = (CPU_Rank == CPU_Number-1 ? 1 : 0);
  setfpe ();			/* Control behavior for floating point
				   exceptions trapping */
  if (argc == 1) PrintUsage (argv[0]);
  strcpy (ParameterFile, "");
  for (i = 1; i < argc; i++) {
    if (*(argv[i]) == '-') {
      if (strspn (argv[i], "-scndovtprfamzib0123456789") != strlen (argv[i]))
	PrintUsage (argv[0]);
      if (strchr (argv[i], 'n'))
	disable = YES;
      if (strchr (argv[i], 'v'))
	verbose = YES;
      if (strchr (argv[i], 't'))
	TimeInfo = YES;
      if (strchr (argv[i], 'c'))
	SloppyCFL = YES;
      if (strchr (argv[i], 'r'))
	TestRadii = YES;
      if (strchr (argv[i], 'p'))
	Profiling = YES;
      if (strchr (argv[i], 'd'))
	debug = YES;
      if (strchr (argv[i], 'b'))
	CentrifugalBalance = YES;
      if (strchr (argv[i], 'm'))
	Merge = YES;
      if (strchr (argv[i], 'a'))
	MonitorIntegral = YES;
      if (strchr (argv[i], 'z'))
	FakeSequential = YES;
      if (strchr (argv[i], 'i'))
	StoreSigma = YES;
      if (strchr (argv[i], '0'))
	OnlyInit = YES;
      if ((argv[i][1] >= '1') && (argv[i][1] <= '9')) {
	GotoNextOutput = YES;
	StillWriteOneOutput = (int)(argv[i][1]-'0');
      }
      if (strchr (argv[i], 's')) {
	Restart = YES;
	i++;
	NbRestart = atoi(argv[i]);
	if ((NbRestart < 0)) {
	  masterprint ("Incorrect restart number\n");
	  PrintUsage (argv[0]);
	}
      }
      if (strchr (argv[i], 'o')) {
	OverridesOutputdir = YES;
	i++;
	sprintf (NewOutputdir, "%s", argv[i]);
      } else {
	if (strchr (argv[i], 'f')) {
	  i++;
	  ScalingFactor = atof(argv[i]);
	  masterprint ("Scaling factor = %g\n", ScalingFactor);
	  if ((ScalingFactor <= 0)) {
	    masterprint ("Incorrect scaling factor\n");
	    PrintUsage (argv[0]);
	  }
	}
      }
    }
    else strcpy (ParameterFile, argv[i]);
  }
  if ((StoreSigma) && !(Restart)) {
    mastererr ("You cannot use tabulated surface density\n");
    mastererr ("in a non-restart run. Aborted\n");
    prs_exit (0);
  }
  if (ParameterFile[0] == 0) PrintUsage (argv[0]);
  ReadVariables (ParameterFile);
  sys = InitPlanetarySystem (PLANETCONFIG);
  ListPlanets (sys);
  SplitDomain ();
  FillPolar1DArrays2D1D ();  // was in the following if.
  if (TestRadii == YES) {
    prs_exit (0);
  }
  if (verbose == YES) 
    TellEverything ();
  if (disable == YES)
    prs_exit (0);
  DumpSources (argc, argv);
  masterprint ("Allocating arrays...");
  fflush(stdout);
  gas_density        = CreatePolarGrid(NRAD, NSEC, "dens");
  gas_v_rad          = CreatePolarGrid(NRAD, NSEC, "vrad");
  gas_v_theta        = CreatePolarGrid(NRAD, NSEC, "vtheta");
  gas_label          = CreatePolarGrid(NRAD, NSEC, "label");
  gas_density_1D     = CreatePolarGrid1D(NRAD1D, "dens1D");
  gas_v_rad_1D       = CreatePolarGrid1D(NRAD1D, "vrad1D");
  gas_v_theta_1D     = CreatePolarGrid1D(NRAD1D, "vtheta1D");
  gas_label_1D       = CreatePolarGrid1D(NRAD1D, "label1D");
  masterprint (" done.\n");
  fflush(stdout);
  OmegaFrame = OMEGAFRAME;
  LostMass = LostLabel = 0.;
  if (Corotating == YES) OmegaFrame = GetPsysInfo (sys, FREQUENCY);
  Initialization2D1D (gas_density, gas_v_rad, gas_v_theta, gas_label, gas_density_1D, gas_v_rad_1D, gas_v_theta_1D, gas_label_1D);
  InitComputeAccel ();
  if (Restart == YES) {
    begin_i         = NbRestart * NINTERM;
    RestartPlanetarySystem (NbRestart, sys);
    LostMass = GetfromPlanetFile (NbRestart, 7, 0); /* 0 refers to planet #0 = star */
    PhysicalTime  = GetfromPlanetFile (NbRestart, 8, 0);
    OmegaFrame  = GetfromPlanetFile (NbRestart, 9, 0);
  } else {			/* We initialize 'planet[i].dat' file */
    EmptyPlanetSystemFile (sys);
  }
  MultiplyPolarGridbyConstant (gas_density, ScalingFactor);
  MultiplyPolarGridbyConstant1D (gas_density_1D, ScalingFactor);
  if (Restart != YES)
    BarycenterConservation(sys,gas_density,gas_v_rad,gas_v_theta);
  for (i = begin_i; i <= NTOT; i++) {
    InnerOutputCounter++;
    if (InnerOutputCounter == 1) {
      InnerOutputCounter = 0;
      WriteBigPlanetSystemFile (sys, TimeStep);
      UpdateLog (sys, gas_density, TimeStep, PhysicalTime);
    }
    if (NINTERM * (TimeStep = (i / NINTERM)) == i) {	/* Outputs are done here */
      TimeToWrite = YES;
      SendOutput2D1D (TimeStep, gas_density, gas_v_rad, gas_v_theta, gas_label, gas_density_1D, gas_v_rad_1D, gas_v_theta_1D, gas_label_1D);
      WritePlanetSystemFile (sys, TimeStep);
      if ((OnlyInit) || ((GotoNextOutput) && (!StillWriteOneOutput))) {
	MPI_Finalize();
	return 0;
      }
      StillWriteOneOutput--;
      if (TimeInfo == YES)	/* Time monitoring is done here */
	GiveTimeInfo (TimeStep);
    }
    else {
      TimeToWrite = NO;
    }
				/* Algorithm loop begins here */

				/***********************/
				/* Hydrodynamical Part */
				/***********************/
    InitSpecificTime (Profiling, &t_Hydro, "Eulerian Hydro algorithms");
    AlgoGas2D1D (gas_density, gas_v_rad, gas_v_theta, gas_label, gas_density_1D, gas_v_rad_1D, gas_v_theta_1D, gas_label_1D, sys);
    GiveSpecificTime (Profiling, t_Hydro);
    SolveOrbits (sys);
    if (MonitorIntegral == YES) {
      masterprint ("Gas Momentum in 2D-Grid   : %.18g\n", GasMomentum (gas_density, gas_v_theta));
      masterprint ("Gas total Mass in 2D-Grid : %.18g\n", GasTotalMass (gas_density));
      masterprint ("Gas Momentum in 1D-Grid   : %.18g\n", GasMomentum1D(gas_density_1D,gas_v_theta_1D));
      masterprint ("Gas total Mass in 1D-Grid : %.18g\n", GasTotalMass1D (gas_density_1D));
    }
  }
  FreePlanetary (sys);
  MPI_Finalize ();
  return 0;
}
 
