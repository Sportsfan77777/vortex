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
#define  VISCOSITY               SUPER_TIME_STEPPING
#define  ROTATING_FRAME          YES

/* -- user-defined parameters (labels) -- */

#define  P_Mstar                 0
#define  P_Mdisk                 1
#define  P_Mplanet               2
#define  P_Viscosity             3
#define  P_AspectRatio           4
#define  P_DensityPower          5
#define  P_TemperaturePower      6
#define  P_Sigma0                7
#define  P_SmoothingLength       8
#define  P_SmoothingType         9
#define  P_IndirectTerm          10
#define  P_ViscosityType         11
#define  P_BaseViscosity         12
#define  P_MaxViscosity          13
#define  P_ViscRampCenter        14
#define  P_ViscRampWidth         15

/* [Beg] user-defined constants (do not change this line) */

#define  UNIT_LENGTH             (5.2*CONST_au)
#define  UNIT_DENSITY            (CONST_Msun/(UNIT_LENGTH*UNIT_LENGTH*UNIT_LENGTH))
#define  UNIT_VELOCITY           (sqrt(CONST_G*g_inputParam[P_Mstar]*CONST_Msun/UNIT_LENGTH)/(2.*CONST_PI))

/* [End] user-defined constants (do not change this line) */

/* -- supplementary constants (user editable) -- */ 

#define  INITIAL_SMOOTHING   NO
#define  WARNING_MESSAGES    YES
#define  PRINT_TO_FILE       NO
#define  INTERNAL_BOUNDARY   NO
#define  SHOCK_FLATTENING    NO
#define  CHAR_LIMITING       NO
#define  LIMITER             VANLEER_LIM
