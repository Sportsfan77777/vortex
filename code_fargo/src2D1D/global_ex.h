/********************************/
/*                              */
/* This file is created         */
/* automatically during         */
/* compilation. Do not edit.    */
/* See perl script              */
/* "varparser.pl" for details   */
/*                              */
/********************************/
extern int CPU_Rank;
extern int CPU_Number;
extern boolean CPU_Master, IAmTheFirst, IAmTheLast;
extern int IMIN;
extern int IMAX;
extern int Zero_or_active;
extern int Max_or_active;
extern int One_or_active;
extern int MaxMO_or_active;		/* MO: Minus One */
extern int GLOBALNRAD;
extern int NRAD1D, IINNER;
extern real Rinf[MAX1D], Rsup[MAX1D], Rmed[MAX1D], Surf[MAX1D];
extern real InvRmed[MAX1D], InvSurf[MAX1D], InvDiffRmed[MAX1D];
extern real InvDiffRsup[MAX1D], InvRinf[MAX1D], Radii[MAX1D], GlobalRmed[MAX1D];
extern real Rinf1D[MAX1D], Rsup1D[MAX1D], Rmed1D[MAX1D], Surf1D[MAX1D];
extern real InvRmed1D[MAX1D], InvSurf1D[MAX1D], InvDiffRmed1D[MAX1D];
extern real InvDiffRsup1D[MAX1D], InvRinf1D[MAX1D], Radii1D[MAX1D], GlobalRmed1D[MAX1D];
extern real SigmaMed1D[MAX1D], SigmaInf1D[MAX1D], MassTaper;
extern real SigmaMed[MAX1D], SigmaInf[MAX1D];
extern real Cosinus[MAX1D], Sinus[MAX1D];
extern real OmegaFrame, PhysicalTime;
extern int TimeStep;
extern real SOUNDSPEED[MAX1D];
extern real SOUNDSPEED1D[MAX1D];
extern real GLOBAL_SOUNDSPEED[MAX1D];
extern boolean Merge, AdvecteLabel, FakeSequential, MonitorIntegral, debug, OnlyInit;
extern boolean	GotoNextOutput, StoreSigma, ViscosityAlpha, RocheSmoothing;
extern boolean CentrifugalBalance, ExcludeHill, SloppyCFL;
extern MPI_Status stat;
extern PolarGrid *CellAbscissa, *CellOrdinate;
extern PolarGrid *RhoStar, *RhoInt, *Potential;
extern PolarGrid1D *RhoStar1D, *RhoInt1D, *Potential1D;
extern PolarGrid *Potplanet;
extern boolean LogGrid;
extern boolean OverridesOutputdir;
extern char NewOutputdir[1024];
