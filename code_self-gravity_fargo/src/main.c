#include "mp.h"

boolean         TimeToWrite, Restart = NO, OpenInner = NO;
int             TimeStep = 0, begin_i = 0, NbRestart = 0, verbose = NO;
int             dimfxy = 11;
static int      InnerOutputCounter=0, StillWriteOneOutput;
extern real     LostMass;
extern boolean  Corotating;
extern boolean  SelfGravity, SGZeroMode, Adiabatic;
real            ScalingFactor = 1.0;

int
main(argc, argv)
     int argc;
     char *argv[];
{
  PolarGrid   *gas_density;
  PolarGrid   *gas_v_rad; 
  PolarGrid   *gas_v_theta; 
  PolarGrid   *gas_energy; 
  PolarGrid   *gas_label;
  int          i;
  real         foostep = 0.;
  boolean      disable = NO, TimeInfo = NO, Profiling = NO;
  boolean      Stockholm = NO, updatevelocities = NO;
  TimeProcess  t_Hydro;
  char         ParameterFile[256];
  PlanetarySystem *sys;
  Force *force;
  
  MPI_Init (&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &CPU_Rank);
  MPI_Comm_size (MPI_COMM_WORLD, &CPU_Number);
  CPU_Master = (CPU_Rank == 0 ? 1 : 0);
  setfpe ();  /* Control behavior for floating point
		 exceptions trapping (default is not to do anything) */
  if (argc == 1) PrintUsage (argv[0]);
  strcpy (ParameterFile, "");
  for (i = 1; i < argc; i++) {
    if (*(argv[i]) == '-') {
      if (strspn (argv[i], "-secndovtpfamzib0123456789") != strlen (argv[i]))
	PrintUsage (argv[0]);
      if (strchr (argv[i], 'n'))
	disable = YES;
      if (strchr (argv[i], 'e'))
	Stockholm = YES;
      if (strchr (argv[i], 'v'))
	verbose = YES;
      if (strchr (argv[i], 't'))
	TimeInfo = YES;
      if (strchr (argv[i], 'c'))
	SloppyCFL = YES;
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
      if (strchr (argv[i], 'i')) {
	StoreSigma = YES;
	if (Adiabatic)
	  StoreEnergy = YES;
      }
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
  if ( (StoreSigma || StoreEnergy) && !(Restart)) {
    mastererr ("You cannot use tabulated surface density\n");
    mastererr ("or surface internal energy in a non-restart run.\n");
    mastererr ("Aborted\n");
    prs_exit (0);
  }
  if (ParameterFile[0] == 0) PrintUsage (argv[0]);
  ReadVariables (ParameterFile);
  SplitDomain ();
  if (verbose == YES) 
    TellEverything ();
  if (disable == YES)
    prs_exit (0);
  DumpSources (argc, argv);
  masterprint ("Allocating arrays...");
  fflush (stdout);
  gas_density        = CreatePolarGrid(NRAD, NSEC, "dens");
  gas_v_rad          = CreatePolarGrid(NRAD, NSEC, "vrad");
  gas_v_theta        = CreatePolarGrid(NRAD, NSEC, "vtheta");
  gas_energy         = CreatePolarGrid(NRAD, NSEC, "energy");
  gas_label          = CreatePolarGrid(NRAD, NSEC, "label");
  masterprint ("done.\n");
  FillPolar1DArrays ();
  force = AllocateForce (dimfxy);
  
  /* Here planets are initialized feeling star potential but they do
     not feel disk potential */
  sys = InitPlanetarySystem (PLANETCONFIG);
  
  /* Gas density initialization */
  InitGasDensity (gas_density);
    
  /* If energy equation is taken into account, we initialize the gas
     thermal energy */
  if ( Adiabatic )
    InitGasEnergy (gas_energy);

  if ( SelfGravity ) {
    /* If SelfGravity = YES or Z, planets are initialized feeling disk
       potential. Only the surface density is required to calculate
       the radial self-gravity acceleration. The disk radial and
       azimutal velocities are not updated */
    compute_selfgravity (gas_density, gas_v_rad, gas_v_theta, foostep, updatevelocities);
    init_planetarysys_withSG (sys);
  }
  ListPlanets (sys);
  OmegaFrame = OMEGAFRAME;
  if (Corotating == YES) OmegaFrame = GetPsysInfo (sys, FREQUENCY);
  
  /* Only gas velocities remain to be initialized */
  Initialization (gas_density, gas_v_rad, gas_v_theta, gas_energy, gas_label, sys);

  /* Initial gas_density is used to compute the circumplanetary mass
     with initial density field */
  mdcp0 = CircumPlanetaryMass (gas_density, sys);
  
  if (Restart == YES) {
    begin_i         = NbRestart * NINTERM;
    RestartPlanetarySystem (NbRestart, sys);
    LostMass = GetfromPlanetFile (NbRestart, 8, 0); /* 0 refers to planet #0 */
    PhysicalTime  = GetfromPlanetFile (NbRestart, 9, 0);
    OmegaFrame  = GetfromPlanetFile (NbRestart, 10, 0);
  } else {			/* We initialize 'planet[i].dat' file */
    EmptyPlanetSystemFile (sys);
  }
  if (MonitorIntegral == YES)
    CheckMomentumConservation (gas_density, gas_v_theta, sys);
  PhysicalTimeInitial = PhysicalTime;
  MultiplyPolarGridbyConstant (gas_density, ScalingFactor);
  for (i = begin_i; i <= NTOT; i++) {
    InnerOutputCounter++;
    if (InnerOutputCounter == 1) {
      InnerOutputCounter = 0;
      WriteBigPlanetSystemFile (sys, TimeStep);
      UpdateLog (force, sys, gas_density, gas_energy, TimeStep, PhysicalTime, dimfxy);
      if (Stockholm == YES)
	UpdateLogStockholm (sys, gas_density, gas_energy, TimeStep, PhysicalTime);
    }
    if (NINTERM * (TimeStep = (i / NINTERM)) == i) {
      /* Outputs are done here */
      TimeToWrite = YES;
      SendOutput (TimeStep, gas_density, gas_v_rad, gas_v_theta, gas_energy, gas_label);
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
    AlgoGas (force, gas_density, gas_v_rad, gas_v_theta, gas_energy, gas_label, sys);
    GiveSpecificTime (Profiling, t_Hydro);
    SolveOrbits (sys);
    
    if (MonitorIntegral == YES) {
      CheckMomentumConservation (gas_density, gas_v_theta, sys);
      masterprint ("Gas Momentum   : %.18g\n", GasMomentum (gas_density, gas_v_theta));
      masterprint ("Gas total Mass : %.18g\n", GasTotalMass (gas_density));
      masterprint ("Gas total Energy : %.18g\n", GasTotalEnergy (gas_density, gas_v_rad, gas_v_theta, gas_energy));
    }
  }
  FreePlanetary (sys);
  FreeForce (force);
  if ( SelfGravity && !SGZeroMode ) {
    rfftwnd_mpi_destroy_plan(SGP_fftplan_forward);
    rfftwnd_mpi_destroy_plan(SGP_fftplan_backward);
  }
  MPI_Finalize ();
  return 0;
}
