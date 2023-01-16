/* 

   Note for the user: 

   This file contains the definition of the collision matrix used to
   solve the drag force implictly. This file is added as an #include
   in both the collisions.c and collisions_template.cu files.

*/

// first index = cols, second index = rows,
// e.g. m[p+o*NFLUIDS] => p col and o row.

idm = lxm*id1 + lym*id2 + lzm*id3;

#if defined(STOKESNUMBER) || defined(FIXPARTICLESIZE)

#ifdef CYLINDRICAL
bigr = ymed(j)*(id1+id3) + ymin(j)*id2;
#endif

#ifdef SPHERICAL
bigr = ymed(j)*(sin(zmed(k)*id1 + zmin(k)*id3)) + ymin(j)*sin(zmed(k))*id2; 
#endif

omega = sqrt(G*MSTAR/bigr/bigr/bigr);
h = ASPECTRATIO*pow(bigr/R0, FLARINGINDEX)*bigr;

#ifdef ISOTHERMAL
//cs = omega*h;
cs = 0.5*(e[0][l] + e[0][idm]);
#endif

#ifdef ADIABATIC
cs = 0.5*(sqrt((GAMMA-1.0)*e[0][l]/rho[0][l]) + sqrt((GAMMA-1.0)*e[0][idm]/rho[0][idm]))
#endif

#endif

#ifdef CONSTANTSTOKESNUMBER
#ifdef SHEARINGBOX
omega = OMEGAFRAME;
#else
omega = 1.0;
#endif
#endif

//Enforce minimum dust-to-gas ratio. Temporary fix here.
for(o=1; o<NFLUIDS; o++){
  rho[o][l] = rho[o][l]/rho[0][l] > 1.0e-16 ? rho[o][l]:1.0e-16*rho[0][l];
 }

// In the implementation below, alpha --> 1/St

for (o=0; o<NFLUIDS; o++) {
  for (p=0; p<NFLUIDS; p++) {
    
    rho_p  = 0.5*(rho[p][l] + rho[p][idm]);
    
    // off-diagonal elements
    if (p != o) {
      
      rho_o  = 0.5*(rho[o][l] + rho[o][idm]);
      
      /* In the line below, the collision term should be
	 alpha[o+p*NFLUIDS], however, we use alpha[p+o*NFLUIDS] to
	 have the possibility of disabling feedback if necessary.*/      
      
#if defined(STOKESNUMBER) || defined(CONSTANTSTOKESNUMBER)
      if ( p > o )  m[p+o*NFLUIDS] = -dt*omega*alpha[p+o*NFLUIDS]*rho_p/rho_o;
      else          m[p+o*NFLUIDS] = -dt*omega*alpha[p+o*NFLUIDS];
#endif
#ifdef CONSTANTDRAG
      m[p+o*NFLUIDS] = -dt*alpha[p+o*NFLUIDS]/rho_o;
#endif

#ifdef FIXPARTICLESIZE
     m[p+o*NFLUIDS] = -dt*alpha[p+o*NFLUIDS]*cs*rho_p;
#endif
    }
    
    // diagonal elements
    else {
      
      /* We now compute the sum present in the diagonal elements.
	 (**) The sum is factorized by dt*omega  */
      
      sum = 0.0;
      for (q=0; q<NFLUIDS; q++) {
	
	//Element pp not included
	if (q != p){
	  
	  rho_q  = 0.5*(rho[q][l] + rho[q][idm]);
	  
#if defined(STOKESNUMBER) || defined(CONSTANTSTOKESNUMBER)
	  
	  /* In the line below, the collision term should be
	     alpha[p+q*NFLUIDS], however, we use alpha[q+p*NFLUIDS] to
	     have the possibility of disabling feedback if necessary.*/
	  
	  if( q > p ) sum += alpha[q+p*NFLUIDS]*rho_q/rho_p;
	  else        sum += alpha[q+p*NFLUIDS];
#endif
#ifdef CONSTANTDRAG
	  sum += alpha[q+p*NFLUIDS];
#endif	  
#ifdef FIXPARTICLESIZE
          sum += alpha[q+p*NFLUIDS]*rho_q;
#endif
	}
      }
      
#if defined(STOKESNUMBER) || defined(CONSTANTSTOKESNUMBER)
      m[p+p*NFLUIDS] = 1.0 + dt*omega*sum; //The factors were not present in the sum (see **)
#endif
      
#ifdef CONSTANTDRAG
      m[p+p*NFLUIDS] = 1.0 + dt*sum/rho_p; //The factors were not present in the sum (see **)
#endif

#ifdef FIXPARTICLESIZE
      m[p+p*NFLUIDS] = 1.0 + dt*cs*sum;
#endif
    }
  }
  b[o] = velocities_input[o][l];
 }
