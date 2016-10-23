#define  PHYSICS                 HD
#define  DIMENSIONS              3
#define  COMPONENTS              3
#define  GEOMETRY                SPHERICAL
#define  BODY_FORCE              POTENTIAL
#define  COOLING                 NO
#define  INTERPOLATION           LINEAR
#define  TIME_STEPPING           RK2
#define  DIMENSIONAL_SPLITTING   NO
#define  NTRACER                 0
#define  USER_DEF_PARAMETERS     28

/* -- physics dependent declarations -- */

#define    EOS                     ISOTHERMAL
#define    ENTROPY_SWITCH          NO
#define    THERMAL_CONDUCTION      NO
#define    VISCOSITY               SUPER_TIME_STEPPING
#define    ROTATING_FRAME          NO

/* -- pointers to user-def parameters -- */

#define  smallq             0
#define  smallp             1
#define  smallh             2
#define  r0                 3
#define  rout               4
#define  gmma               5
#define  mplanet            6
#define  amp                7
#define  bump_width         8
#define  qout               9
#define  rmin               10
#define  pert_amp           11
#define  damp_in            12
#define  damp_out           13
#define  tdamp              14
#define  soft               15
#define  planet_on          16
#define  switch_on          17
#define  nu                 18
#define  nu_alpha           19
#define  visc_jump_amp      20
#define  visc_jump_height   21
#define  visc_jump_width    22
#define  visc_jump_psimax   23
#define  visc_jump_H        24
#define  rad_obc            25
#define  damp_up            26
#define  visc_rtrans        27

/* -- supplementary constants (user editable) -- */ 

#define  INITIAL_SMOOTHING     NO
#define  WARNING_MESSAGES      YES
#define  PRINT_TO_FILE         NO
#define  INTERNAL_BOUNDARY     YES
#define  SHOCK_FLATTENING      NO
#define  ARTIFICIAL_VISCOSITY  NO
#define  CHAR_LIMITING         NO
#define  LIMITER               VANLEER_LIM
#define  STS_nu                0.01
