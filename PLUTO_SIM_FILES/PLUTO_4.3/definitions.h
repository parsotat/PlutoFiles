#define  PHYSICS                        RHD
#define  DIMENSIONS                     2
#define  COMPONENTS                     2
#define  GEOMETRY                       SPHERICAL
#define  BODY_FORCE                     NO
#define  FORCED_TURB                    NO
#define  COOLING                        NO
#define  RECONSTRUCTION                 PARABOLIC
#define  TIME_STEPPING                  RK2
#define  DIMENSIONAL_SPLITTING          NO
#define  NTRACER                        0
#define  USER_DEF_PARAMETERS            14

/* -- physics dependent declarations -- */

#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 NO

/* -- user-defined parameters (labels) -- */

#define  JET_ON                         0
#define  PRESS_JET                      1
#define  RHO_JET                        2
#define  VEL_JET                        3
#define  THETA_JET                      4
#define  R_MIN                          5
#define  R_STAR                         6
#define  PRESS_AMBIENT                  7
#define  RHO_AMBIENT                    8
#define  ADIABATIC_INDEX                9
#define  JET_PERIOD                     10
#define  JET_FRAC_ON                    11
#define  JET_LUMI_INIT_FRAC             12
#define  JET_LUMI_FRAC_DEC              13


/* [Beg] user-defined constants (do not change this line) */

#define  CHOMBO_LOGR                    YES
#define  RECONSTRUCT_4VEL               YES
#define  SHOCK_FLATTENING               MULTID
#define  CHAR_LIMITING                  YES
#define  EPS_PSHOCK_FLATTEN             2.5
#define  ARTIFICIAL_VISC                0.1
#define  UNIT_LENGTH                    1e9
#define  UNIT_VELOCITY                  CONST_c
#define  UNIT_DENSITY                   1.0
#define  CHOMBO_REF_VAR                 RHO
#define  INTERNAL_BOUNDARY              YES
/* [End] user-defined constants (do not change this line) */
