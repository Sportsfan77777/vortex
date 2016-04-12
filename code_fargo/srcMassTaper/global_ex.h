/** \file global_ex.h

This file is created      
automatically during      
compilation from global.h. Do not edit. 
See perl script           
"varparser.pl" for details

\file global.h

Declares all global variables.
Used to construct automatically
the file global_ex.h. The file
global.h cannot contain any comment,
as it would not be parsed correctly
by varparser.pl
*/                          

extern int CPU_Rank;
extern int CPU_Number;
extern boolean CPU_Master;
extern int IMIN;
extern int IMAX;
extern int Zero_or_active;
extern int Max_or_active;
extern int One_or_active;
extern int MaxMO_or_active;		/* MO: Minus One */
extern int GLOBALNRAD;
extern real Rinf[MAX1D], Rsup[MAX1D], Rmed[MAX1D], Surf[MAX1D];
extern real InvRmed[MAX1D], InvSurf[MAX1D], InvDiffRmed[MAX1D];
extern real InvDiffRsup[MAX1D], InvRinf[MAX1D], Radii[MAX1D], GlobalRmed[MAX1D];
extern real SigmaMed[MAX1D], SigmaInf[MAX1D], MassTaper;
extern real OmegaFrame, PhysicalTime, PhysicalTimeInitial;
extern int TimeStep;
extern real SOUNDSPEED[MAX1D];
extern real GLOBAL_SOUNDSPEED[MAX1D];
extern boolean Merge, AdvecteLabel, FakeSequential, MonitorIntegral, debug, OnlyInit;
extern boolean	GotoNextOutput, StoreSigma, ViscosityAlpha, RocheSmoothing;
extern boolean CentrifugalBalance, ExcludeHill, SloppyCFL, Stockholm;
extern MPI_Status fargostat;
extern PolarGrid *CellAbscissa, *CellOrdinate;
extern PolarGrid *RhoStar, *RhoInt, *Potential;
extern boolean LogGrid;
extern boolean OverridesOutputdir;
extern char NewOutputdir[1024];
extern boolean FakeAccretion;
