#define  PHYSICS                 HD
#define  DIMENSIONS              2
#define  COMPONENTS              2
#define  GEOMETRY                POLAR
#define  BODY_FORCE              POTENTIAL
#define  COOLING                 NO
#define  INTERPOLATION           LINEAR
#define  TIME_STEPPING           RK2
#define  DIMENSIONAL_SPLITTING   NO
#define  NTRACER                 0
#define  USER_DEF_PARAMETERS     6
#define  USER_DEF_CONSTANTS      3

/* -- physics dependent declarations -- */

#define  EOS                     ISOTHERMAL
#define  ENTROPY_SWITCH          NO
#define  THERMAL_CONDUCTION      NO
#define  VISCOSITY               SUPER_TIME_STEPPING
#define  ROTATING_FRAME          YES

/* -- user-defined parameters (labels) -- */

#define  Mstar                   0
#define  Mdisk                   1
#define  Mplanet                 2
#define  Viscosity               3
#define  DensityPower            4
#define  SoundSpeedPower         5

/* -- user-defined symbolic constants -- */

#define  UNIT_LENGTH             (1.0*CONST_au)
#define  UNIT_DENSITY            (CONST_Msun/(UNIT_LENGTH*UNIT_LENGTH*UNIT_LENGTH))
#define  UNIT_VELOCITY           sqrt(CONST_G*g_inputParam[Mstar]*CONST_Msun/UNIT_LENGTH)/(2.*CONST_PI)

/* -- supplementary constants (user editable) -- */ 

#define  INITIAL_SMOOTHING      NO
#define  WARNING_MESSAGES       YES
#define  PRINT_TO_FILE          NO
#define  INTERNAL_BOUNDARY      NO
#define  SHOCK_FLATTENING       NO
#define  ARTIFICIAL_VISCOSITY   NO
#define  CHAR_LIMITING          NO
#define  LIMITER                VANLEER_LIM
#define  STS_nu                 0.01
