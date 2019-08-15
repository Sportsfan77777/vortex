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
extern boolean CPU_Master;
extern int CPU_Next, CPU_Prev, CPU_Highest;
extern int CPU_Friend, CPU_NoFriend;
extern real *dens_friend;
extern real *SGP_buffft_Accr_friend, *SGP_buffft_Acct_friend;
extern real *ffttohydro_transfer, *ffttohydro_transfer_friend;
extern int local_Nx, local_i_start, local_i_start_friend, total_local_size_friend, local_Nx_friend;
extern int local_Ny_after_transpose, local_j_start_after_transpose;
extern int total_local_size, ifront, Zero_or_active_friend;
extern int transfer_size, transfer_size_friend;
extern int hydro_totalsize, active_hydro_totalsize, active_hydro_totalsize_friend;  
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
extern real EnergyMed[MAX1D], QplusMed[MAX1D], CoolingTimeMed[MAX1D];
extern real OmegaFrame, PhysicalTime, PhysicalTimeInitial;
extern int TimeStep;
extern real HillRadius, mdcp, mdcp0, exces_mdcp;
extern real GLOBAL_bufarray[MAX1D];
extern boolean Merge, AdvecteLabel, FakeSequential, MonitorIntegral, debug, OnlyInit;
extern boolean	GotoNextOutput, StoreSigma, StoreEnergy, ViscosityAlpha, RocheSmoothing;
extern boolean CentrifugalBalance, ExcludeHill, SloppyCFL;
extern MPI_Status stat;
extern PolarGrid *CellAbscissa, *CellOrdinate;
extern PolarGrid *RhoStar, *RhoInt, *Potential, *Pressure, *SoundSpeed, *Temperature;
extern PolarGrid *DivergenceVelocity, *TAURR, *TAUPP, *TAURP, *Qplus;
extern boolean LogGrid;
extern boolean OverridesOutputdir;
extern char NewOutputdir[1024];
extern boolean FakeAccretion;
