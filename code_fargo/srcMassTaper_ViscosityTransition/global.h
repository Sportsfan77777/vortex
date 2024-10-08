int CPU_Rank;
int CPU_Number;
boolean CPU_Master;
int IMIN;
int IMAX;
int Zero_or_active;
int Max_or_active;
int One_or_active;
int MaxMO_or_active;		/* MO: Minus One */
int GLOBALNRAD;
real Rinf[MAX1D], Rsup[MAX1D], Rmed[MAX1D], Surf[MAX1D];
real InvRmed[MAX1D], InvSurf[MAX1D], InvDiffRmed[MAX1D];
real InvDiffRsup[MAX1D], InvRinf[MAX1D], Radii[MAX1D], GlobalRmed[MAX1D];
real SigmaMed[MAX1D], SigmaInf[MAX1D], MassTaper;
real OmegaFrame, PhysicalTime=0.0, PhysicalTimeInitial;
int TimeStep=0;
real SOUNDSPEED[MAX1D];
real GLOBAL_SOUNDSPEED[MAX1D];
boolean Merge, AdvecteLabel, FakeSequential, MonitorIntegral, debug, OnlyInit;
boolean	GotoNextOutput, StoreSigma, ViscosityAlpha, RocheSmoothing;
boolean CentrifugalBalance, ExcludeHill, SloppyCFL, Stockholm, DeadZone, FullDeadZone, EvolvingDeadZone;
MPI_Status fargostat;
PolarGrid *CellAbscissa, *CellOrdinate;
PolarGrid *RhoStar, *RhoInt, *Potential;
boolean LogGrid;
boolean OverridesOutputdir;
char NewOutputdir[1024];
boolean FakeAccretion;
boolean SinSquaredTaper, ParabolaTaper;
boolean IncludeAccretionAngularMomentum;
