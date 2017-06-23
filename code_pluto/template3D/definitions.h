#define  PHYSICS                 HD
#define  DIMENSIONS              2
#define  COMPONENTS              3
#define  GEOMETRY                SPHERICAL
#define  BODY_FORCE              POTENTIAL
#define  COOLING                 NO
#define  RECONSTRUCTION          LINEAR
#define  TIME_STEPPING           RK2
#define  DIMENSIONAL_SPLITTING   NO
#define  NTRACER                 0
#define  USER_DEF_PARAMETERS     16

/* -- physics dependent declarations -- */

#define  EOS                     ISOTHERMAL
#define  ENTROPY_SWITCH          NO
#define  THERMAL_CONDUCTION      NO
#define  VISCOSITY               YES
#define  ROTATING_FRAME          YES

/* -- user-defined parameters (labels) -- */

#define  P_Mstar                 0
#define  P_Mplanet               1
#define  P_Viscosity             2
#define  P_AspectRatio           3
#define  P_DensityPower          4
#define  P_TemperaturePower      5
#define  P_Sigma0                6
#define  P_SmoothingLength       7
#define  P_SMOOTH_SCALE_HEIGHT   8
#define  P_SMOOTH_HILL_RADIUS    9
#define  P_INDIRECT_TERM        10
#define  P_ViscosityType        11
#define  P_BaseViscosity        12
#define  P_MaxViscosity         13
#define  P_ViscosityRampCenter  14
#define  P_ViscosityRampWidth   15

/* [Beg] user-defined constants (do not change this line) */

#define  UNIT_LENGTH             (5.2*CONST_au)
#define  UNIT_DENSITY            (CONST_Msun/(UNIT_LENGTH*UNIT_LENGTH*UNIT_LENGTH))
#define  UNIT_VELOCITY           (sqrt(CONST_G*g_inputParam[P_Mstar]*CONST_Msun/UNIT_LENGTH))

/* [End] user-defined constants (do not change this line) */

/* -- supplementary constants (user editable) -- */ 

#define  INITIAL_SMOOTHING   NO
#define  WARNING_MESSAGES    YES
#define  PRINT_TO_FILE       NO
#define  INTERNAL_BOUNDARY   YES
#define  SHOCK_FLATTENING    NO
#define  CHAR_LIMITING       NO
#define  LIMITER             VANLEER_LIM
