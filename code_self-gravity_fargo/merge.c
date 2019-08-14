#include "mp.h"

void merge (nb)
     int nb;
{
  extern boolean  Write_Density, Write_Velocity, Write_Energy, IsDisk;
  extern boolean  Write_Temperature, Write_DivV, Write_Qplus;
  extern boolean SelfGravity, SGZeroMode;
  boolean bool = NO;
  int i, j, l;
  int one_if_odd;
  char radix[512];
  char command[1024];
  if (!CPU_Master) return;
  message ("Merging output files...");
  for (j = 0; j < 7+(AdvecteLabel == YES); j++) {
    switch (j) {
    case 0:
      strcpy (radix, "dens");
      bool = Write_Density;
      break;
    case 1: 
      strcpy (radix, "vrad");
      bool = Write_Velocity;
      break;
    case 2: 
      strcpy (radix, "vtheta");
      bool = Write_Velocity;
      break;
    case 3: 
      strcpy (radix, "energy");
      bool = Write_Energy;
      break;
    case 4: 
      strcpy (radix, "Temperature");
      bool = Write_Temperature;
      break;
    case 5: 
      strcpy (radix, "DivV");
      bool = Write_DivV;
      break;
    case 6: 
      strcpy (radix, "Qplus");
      bool = Write_Qplus;
      break;
    case 7: 
      strcpy (radix, "label");
      bool = YES;
      break;
    }
    one_if_odd = (CPU_Number%2 == 0 ? 0 : 1);
    if (bool == YES) {
      if ( SelfGravity && !SGZeroMode ) {
	for (i = 0; i < (CPU_Number+one_if_odd)/2; i++) {
	  if ( i != 0 ) {
	    sprintf (command, "cd %s; cat gas%s%d.dat.%05d >> gas%s%d.dat", \
		     OUTPUTDIR, radix, nb, i, radix, nb);
	    system (command);
	  }
	  l = (i + (CPU_Number + one_if_odd)/2)%(CPU_Number + one_if_odd);
	  if ( i != CPU_Highest ) {
	    sprintf (command, "cd %s; cat gas%s%d.dat.%05d >> gas%s%d.dat", \
		     OUTPUTDIR, radix, nb, l, radix, nb);
	    system (command);
	  }
	}
	sprintf (command, "cd %s; rm -f gas%s%d.dat.0*",	\
		 OUTPUTDIR, radix, nb);
	system (command);
      }
      else {
	for (i = 1; i < CPU_Number; i++) {
	  sprintf (command, "cd %s; cat gas%s%d.dat.%05d >> gas%s%d.dat", \
		   OUTPUTDIR, radix, nb, i, radix, nb);
	  system (command);
	}
	sprintf (command, "cd %s; rm -f gas%s%d.dat.0*",	\
		 OUTPUTDIR, radix, nb);
	system (command);
      }
    }
  }
  message ("done\n");
  fflush (stdout);
}
