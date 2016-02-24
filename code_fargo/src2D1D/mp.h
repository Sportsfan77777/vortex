#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include "fondam.h"
#include "types.h"
#include "proto.h"
#ifdef _PARALLEL
#include <mpi.h>
#else
#include "mpi_dummy.h"
#endif
#ifndef __LOCAL
#include "param.h"
#include "global_ex.h"
#else
#include "param_noex.h"
#include "global.h"
#endif
#include <sys/times.h>
#include <unistd.h>
#include <string.h>
#ifndef __APPLE__
#include <malloc.h>
#endif
#ifdef _TRAP_FPE
#include <signal.h>
#include <fenv.h>
#endif
