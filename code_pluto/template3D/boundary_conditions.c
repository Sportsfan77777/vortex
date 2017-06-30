/* ///////////////////////////////////////////////////////////////////// */
/* 
Different types of boundary conditions (This file does nothing.)
*/
/* ///////////////////////////////////////////////////////////////////// */


/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/* 
 *
 *********************************************************************** */
{
  int   i, j, k, nv;
  double *x1, *x2, *x3, R, th, z, OmegaK, v[256];
  
  x1 = grid[IDIR].x;
  x2 = grid[JDIR].x;
  x3 = grid[KDIR].x;

  /// Inner Boundary ///

  if (side == X1_BEG){
    X1_BEG_LOOP(k,j,i){
      #if GEOMETRY == SPHERICAL
         R = x1[i]*sin(x2[j]);
         th = x2[j];
         z = x1[i]*cos(x2[j]);
      #endif
      d->Vc[RHO][k][j][i]   = density3D(R, z);
      d->Vc[VX1][k][j][i]   = 0; radialVelocity_rComponent(R, th, z);
      #if GEOMETRY == CYLINDRICAL
          d->Vc[VX2][k][j][i]   = 0.0; // vz or vtheta
      #elif GEOMETRY == SPHERICAL
          d->Vc[VX2][k][j][i]   = 0; radialVelocity_thetaComponent(R, th, z); // vtheta
      #elif GEOMETRY == POLAR && DIMENSIONS == 3
          d->Vc[VX3][k][j][i]   = 0.0; // vz
      #endif
      d->Vc[iVPHI][k][j][i] = azimuthalVelocity3D(R, z);
      #if EOS == IDEAL
          d->Vc[PRS][k][j][i]  = pressure(R, z);
      #endif
    }
  }

  /// Outer Boundary ///

  if (side == X1_END){
    X1_END_LOOP(k,j,i){
      #if GEOMETRY == SPHERICAL
         R = x1[i]*sin(x2[j]);
         th = x2[j];
         z = x1[i]*cos(x2[j]);
      #endif
      d->Vc[RHO][k][j][i]   = density3D(R, z);
      d->Vc[VX1][k][j][i]   = 0; radialVelocity_rComponent(R, th, z);
      #if GEOMETRY == CYLINDRICAL
          d->Vc[VX2][k][j][i]   = 0.0; // vphi
      #elif GEOMETRY == SPHERICAL
          d->Vc[VX2][k][j][i]   = 0; radialVelocity_thetaComponent(R, th, z); // vtheta
      #elif GEOMETRY == POLAR && DIMENSIONS == 3
          d->Vc[VX3][k][j][i]   = 0.0; // vtheta
      #endif
      d->Vc[iVPHI][k][j][i] = azimuthalVelocity3D(R, z);
      #if EOS == IDEAL
          d->Vc[PRS][k][j][i]  = pressure(R, z);
      #endif
    }
  }

  /// Upper Boundary (z != 0) ///

  if (side == X2_BEG){
    X2_BEG_LOOP(k,j,i){
      #if GEOMETRY == SPHERICAL
         R = x1[i]*sin(x2[j]);
         th = x2[j];
         z = x1[i]*cos(x2[j]);
      #endif
      d->Vc[RHO][k][j][i]   = density3D(R, z);
      d->Vc[VX1][k][j][i]   = 0; radialVelocity_rComponent(R, th, z);
      #if GEOMETRY == CYLINDRICAL
          d->Vc[VX2][k][j][i]   = 0.0; // vz or vtheta
      #elif GEOMETRY == SPHERICAL
          d->Vc[VX2][k][j][i]   = 0; radialVelocity_thetaComponent(R, th, z); // vtheta
      #elif GEOMETRY == POLAR && DIMENSIONS == 3
          d->Vc[VX3][k][j][i]   = 0.0; // vz
      #endif
      d->Vc[iVPHI][k][j][i] = azimuthalVelocity3D(R, z);
      #if EOS == IDEAL
          d->Vc[PRS][k][j][i]  = pressure(R, z);
      #endif
    }
  }

  /// Internal Boundary (External Torque) ///

  if (side == 0) {
    DOM_LOOP(k, j, i) {
      #if GEOMETRY == SPHERICAL
         R = x1[i]*sin(x2[j]);
         th = x2[j];
         z = x1[i]*cos(x2[j]);
      #endif
      d->Vc[VX3][k][j][i] += externalTorque(R, z) * g_dt;
    }
  }

/* original below ************************

  if (side == X1_BEG){
    X1_BEG_LOOP(k,j,i){
      NVAR_LOOP(nv) d->Vc[nv][k][j][i] = d->Vc[nv][k][j][2*IBEG - i - 1];
      d->Vc[VX1][k][j][i] *= -1.0;
      #if GEOMETRY == POLAR
       R = x1[i];
      #elif GEOMETRY == SPHERICAL
       R = x1[i]*sin(x2[j]);
      #endif
      OmegaK = 2.0*CONST_PI/(R*sqrt(R));
      d->Vc[iVPHI][k][j][i] = R*(OmegaK - g_OmegaZ);
#if DUST == YES      
//      NDUST_LOOP(nv) d->Vc[nv][k][j][i] = 0.0;
      d->Vc[VX2_D][k][j][i] = d->Vc[iVPHI][k][j][i];
#endif      
    }
  }

  if (side == X1_END){
    X1_END_LOOP(k,j,i){
      NVAR_LOOP(nv)  d->Vc[nv][k][j][i] = d->Vc[nv][k][j][IEND];
      #if GEOMETRY == POLAR
       R = x1[i];
//       d->Vc[iVR][k][j][i] = 0.0;
      #elif GEOMETRY == SPHERICAL
       R = x1[i]*sin(x2[j]);
       d->Vc[iVR][k][j][i]  = 0.0;
       d->Vc[iVTH][k][j][i] = 0.0;
      #endif
      OmegaK = 2.0*CONST_PI/(R*sqrt(R));
      d->Vc[iVPHI][k][j][i] = R*(OmegaK - g_OmegaZ);
#if DUST == YES      
      d->Vc[VX2_D][k][j][i] = d->Vc[iVPHI][k][j][i];
#endif      
    }
  }
  ********************** */
}