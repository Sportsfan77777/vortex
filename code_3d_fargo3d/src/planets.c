#include "fargo3d.h"

void AccreteOntoPlanets(real dt) {

    //<USER_DEFINED>
    INPUT(Density);
    OUTPUT(Density);
  #ifdef X
    INPUT(Vx);
  #endif
  #ifdef Y
    INPUT(Vy);
  #endif
  #ifdef Z
    INPUT(Vz);
  #endif
  //<\USER_DEFINED>

  int i, j, k, n, i_min, i_max, j_min, j_max, i_f, NbPlanets, count_cells, count1, count2, timestep_a;
  real m, x, y, z, r, r_roche, angle, vx, vy, px, py, distance, dx, dy, MassTaper, surf, core_coef, reduction_factor;
  real negative_taper, negative_time;
  real vrcell, vtcell, vxcell, vycell;
  real facc, facc1, facc2, frac1, frac2;
  real* dens = Density->field_cpu;
  real* vrad = Vy->field_cpu;
  real* vtheta = Vx->field_cpu;

  real dMplanet, dPxPlanet, dPyPlanet, deltaM, temp;
  NbPlanets = Sys->nb;

  real rad_l, az_l;
  char name[256];
  FILE *planet_file;

  count_cells = 0;
  count1 = 0;
  count2 = 0;

  // Return if no accretion
  if (Sys->acc[0] < 0.001) return;

  // Turn printing off
  timestep_a = floor(PhysicalTime / 2.0 / M_PI);
  // sprintf(fn, "accretion%d.dat", timestep);

  // if (access(fn, F_OK) != -1) {
  //   // Already Exists
  //   WriteAccretionFile = NO;
  // }
  // else {
  //   // Create file and open in append mode
  //   if (timestep >= 10) {
  //     file = fopen(fn, "a");
  //     WriteAccretionFile = YES;
  //     printf ("...Accretion file ready...");
  //     fflush (stdout);
  //   }
  // }

  // Calculate mass taper
  if (MASSTAPER == 0.0)
    MassTaper = 1.0;
  else
    MassTaper = (PhysicalTime >= (2.0 * M_PI * MASSTAPER) ? 1.0 : 0.5 * (1.0 - cos(M_PI * PhysicalTime / (2.0 * M_PI * MASSTAPER))));

  if (PhysicalTime < 2.0 * M_PI * NEGATIVESTARTTIME)
    negative_taper = 0.0;
  else {
    negative_time = PhysicalTime - 2.0 * M_PI * NEGATIVESTARTTIME;
    if (negative_time < 0.0)
      negative_taper = 0.0;
    else
      negative_taper = (negative_time >= (2.0 * M_PI * NEGATIVEMASSTAPER) ? 1.0 : .5*(1.0-cos(M_PI*negative_time/(2.0 * M_PI * NEGATIVEMASSTAPER))));
  }

  // Loop through planets
  for (n = 0; n < NbPlanets; n++) {
    dMplanet = dPxPlanet = dPyPlanet = 0.0;
    count1 = count2 = 0;

    // Kley's parameters
    facc = dt * (Sys->acc[n]);
    //facc = (Sys->acc[k]); // Hack-ish test
    facc1 = KLEYFACC1 * facc; // 1.0 / 3.0 * facc;
    facc2 = KLEYFACC2 * facc; // 2.0 / 3.0 * facc;
    frac1 = KLEYFRAC1; // 0.75; (outer region)
    frac2 = KLEYFRAC2; // 0.45; (inner region)

    // Planet parameters
    m = Sys->mass[n] * MassTaper + Sys->accreted_mass[n] - NEGATIVEMASS*negative_taper + 0.000001;
    x = Sys->x[n];
    y = Sys->y[n];
    z = Sys->z[n];

    //reduction_factor = Sys->reduction_factor[n];
    //if (reduction_factor > 1)
       reduction_factor = 1.0; // do not reduce accretion rate
    //if (reduction_factor < 0.1)
    //   reduction_factor = REDUCTION_LIMIT; // anything past this is too much

    r = sqrt(x*x + y*y + z*z);
    r_roche = pow((1.0 / 3.0 * m / MSTAR), (1.0 / 3.0)) * r;
    angle = atan2(y, x);

    vx = Sys->vx[n];
    vy = Sys->vy[n];

    px = m * vx;
    py = m * vy;

    // Indices
    j_min = 0;
    j_max = Ny;
    while ((ymed(j_min) < r - r_roche) && (j_min < Ny)) j_min++; // should be y[j+1]
    while ((ymed(j_max) > r + r_roche) && (j_max > 0)) j_max--; // should be regular y[j]
    j_min -= 2;
    j_max += 2;

    i_min = (int)((real)Nx / 2.0 / M_PI * (angle - 2.0 * r_roche / r)); // does this take into account -pi to pi?
    i_max = (int)((real)Nx / 2.0 / M_PI * (angle + 2.0 * r_roche / r));
    i_min = 0;
    i_max = Nx;

    // Surface density
    surf = 0.0;
    core_coef = 1.0;

    if (COREACCRETION) {
      if (Sys->accreted_mass[n] < Sys->mass[n])
          core_coef = 0.33 + 0.67 * pow(Sys->accreted_mass[n] / Sys->mass[n], 5.0);
    }

    //sprintf (name, "%saccretion%d.dat", OUTPUTDIR, timestep_a);
    //planet_file = fopen_prs (name, "a");

    //fprintf(planet_file, "\nPlanet: <%.4f, %.4f> at %.4f\n", 1000.0 * m, r_roche, r);
    //fprintf(planet_file, "\ni: <%d to %d>; j: <%d to %d> (%.4f to %.4f)\n", i_min, i_max, j_min, j_max, ymed(j_min), ymed(j_max));

    //<MAIN_LOOP>
    #ifdef Z
      for (k=0; k<Nz; k++) {
    #endif
    #ifdef Y
        for (j=j_min; j<j_max; j++) { // radius
    #endif
    #ifdef X
          for (i=i_min; i<i_max; i++) { // azimuth
    #endif
            rad_l = ymed(j);
            az_l = xmed(i);
            //i_f = i;
            //while (i_f <  0)  i_f += Nx;
            //while (i_f >= Nx) i_f -= Nx;

            //xc = abs[l];
            //yc = ord[l];
            dx = x - XC;
            dy = y - YC;

            distance = sqrt(dx*dx + dy*dy);

            count_cells++;
            //fprintf(planet_file, "\nQ %04d: (%.4f, %.4f) (%.4f, %.4f) (%.6f, %.6f, %.6f) ", count_cells, XC, YC, rad_l, az_l, dx, dy, distance);

            vtcell = 0.5 * (vtheta[l] + vtheta[lxp]) + rad_l * OMEGAFRAME;
            vrcell = 0.5 * (vrad[l] + vrad[lyp]);
            vxcell = (vrcell * XC - vtcell * YC) / rad_l;
            vycell = (vrcell * YC + vtcell * XC) / rad_l;

            // Condition 1
            if (distance < frac1 * r_roche) {
              surf = M_PI * (ymed(j+1) * ymed(j+1) - rad_l * rad_l) /(real)Nx; // Set accreted mass here!

              deltaM = facc1 * ACCOFFSET * core_coef * dens[l] * surf * reduction_factor;
              dens[l] *= (1.0 - facc1 * core_coef * reduction_factor);

              dPxPlanet += deltaM * vxcell;
              dPyPlanet += deltaM * vycell;
              dMplanet += deltaM;
              //if ((timestep_a > 62.73 && timestep_a < 62.93) || (timestep_a > 125.55 && timestep_a < 125.75) || (timestep_a > 188.4 && timestep_a < 188.6) || (timestep_a > 251.2 && timestep_a < 251.4))
                //printf("%s", "1");
              //fprintf(planet_file, "\nQ %04d: (%.4f, %.4f) (%.4f, %.4f) (%.6f, %.6f, %.6f) ", count_cells, XC, YC, rad_l, az_l, dx, dy, distance);
                //fprintf(planet_file, "1 ");
              count1++;
              //accreting1++;
            }
            
            // Condition 2
            if (distance < frac2 * r_roche) {
              deltaM = facc2 * ACCOFFSET * core_coef * dens[l] * surf * reduction_factor;
              dens[l] *= (1.0 - facc2 * core_coef * reduction_factor);

              dPxPlanet += deltaM * vxcell;
              dPyPlanet += deltaM * vycell;
              dMplanet += deltaM;
              //if ((timestep_a > 62.73 && timestep_a < 62.93) || (timestep_a > 125.55 && timestep_a < 125.75) || (timestep_a > 188.4 && timestep_a < 188.6) || (timestep_a > 251.2 && timestep_a < 251.4))
                //printf("%s", "2");
                //fprintf(planet_file, "2 ");
              count2++;
              //accreting2++;
            }
            //printf("%d %d", count1, count2);
            //if ((timestep_a > 62.73 && timestep_a < 62.93) || (timestep_a > 125.55 && timestep_a < 125.75) || (timestep_a > 188.4 && timestep_a < 188.6) || (timestep_a > 251.2 && timestep_a < 251.4))
            //  printf("Q");
            // if (WriteAccretionFile && timestep >= 10) {
            //   if (accreting1 || accreting2) {
            //     fprintf (file, "(%.3f, %.3f) - %d %d\n", rad_l, az_l, accreting1, accreting2);
            //   }
            // }
    #ifdef X
          }
    #endif
    #ifdef Y
        }
    #endif
    #ifdef Z
      }
    #endif
    //<\MAIN_LOOP>

    //if ((timestep_a > 62.73 && timestep_a < 62.93) || (timestep_a > 125.55 && timestep_a < 125.75) || (timestep_a > 188.4 && timestep_a < 188.6) || (timestep_a > 251.2 && timestep_a < 251.4))
      //printf("W");
      //fprintf(planet_file, "\nNumber of Accreting Cells: [%d, %d, %d]\n", count_cells, count1, count2);
      //fclose(planet_file);
    //fflush(stdout);

    MPI_Barrier(MPI_COMM_WORLD);

    // They all should reach here at the same time.

    MPI_Allreduce (&dMplanet, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    dMplanet = temp;

    MPI_Allreduce (&dPxPlanet, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    dPxPlanet = temp;
    MPI_Allreduce (&dPyPlanet, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    dPyPlanet = temp;

    //MPI_Allreduce (&count1, &temp, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    //count1 = temp;
    //MPI_Allreduce (&count2, &temp, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    //count2 = temp;

    // if (WriteAccretionFile && timestep >= 10) {
    //   fprintf(file, "Totals: %d %d \n", count1, count2);
    //   fclose(file);
    // }

    px += dPxPlanet;
    py += dPyPlanet;
    //m += dMplanet;

    // Post-processing
    if (FAKEACCRETION) {
          // pass (get rid of disk material, but don't let it affect planet)
      //printf("Fake\n");
    }
    else {
        //Sys->mass[k] = m;
        Sys->accreted_mass[k] += dMplanet;
        //printf("Not Fake %.8f \n", Sys->accreted_mass[k]);
    }

    if (Sys->FeelDisk[n]) {
        //Sys->vx[n] = px / m;
        //Sys->vy[n] = py / m;
    }
  }
}


void ComputeIndirectTerm () {
#ifndef NODEFAULTSTAR
  IndirectTerm.x = -DiskOnPrimaryAcceleration.x;
  IndirectTerm.y = -DiskOnPrimaryAcceleration.y;
  IndirectTerm.z = -DiskOnPrimaryAcceleration.z;
  if (!INDIRECTTERM) {
    IndirectTerm.x = 0.0;
    IndirectTerm.y = 0.0;
    IndirectTerm.z = 0.0;
  }
#else
  IndirectTerm.x = 0.0;
  IndirectTerm.y = 0.0;
  IndirectTerm.z = 0.0;
#endif
}

Force ComputeForce(real x, real y, real z,
		   real rsmoothing, real mass) {
  
  Force Force;

  /* The trick below, which uses VxMed as a 2D temporary array,
     amounts to subtracting the azimuthally averaged density prior to
     the torque evaluation. This has no impact on the torque, but has
     on the angular speed of the planet and is required for a proper
     location of resonances in a non self-gravitating disk. See
     Baruteau & Masset 2008, ApJ, 678, 483 (arXiv:0801.4413) for
     details. */
#ifdef BM08
  ComputeVmed (Total_Density);
  ChangeFrame (-1, Total_Density, VxMed);
#endif
  /* The density is now the perturbed density */
  FARGO_SAFE(_ComputeForce(x, y, z, rsmoothing, mass)); /* Function/Kernel Launcher. */
  /* We restore the total density below by adding back the azimuthal
     average */
#ifdef BM08
  ChangeFrame (+1, Total_Density, VxMed);
#endif

  
#ifdef FLOAT
  MPI_Allreduce (&localforce, &globalforce, 12, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
#else
  MPI_Allreduce (&localforce, &globalforce, 12, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
  
  Force.fx_inner    = globalforce[0];
  Force.fy_inner    = globalforce[1];
  Force.fz_inner    = globalforce[2];
  Force.fx_ex_inner = globalforce[3];
  Force.fy_ex_inner = globalforce[4];
  Force.fz_ex_inner = globalforce[5];
  Force.fx_outer    = globalforce[6];
  Force.fy_outer    = globalforce[7];
  Force.fz_outer    = globalforce[8];
  Force.fx_ex_outer = globalforce[9];
  Force.fy_ex_outer = globalforce[10];
  Force.fz_ex_outer = globalforce[11];

  return Force;

}

Point ComputeAccel(real x, real y, real z,
		   real rsmoothing, real mass) {
  Point acceleration;
  Force force;
  force = ComputeForce (x, y, z, rsmoothing, mass);
  if (EXCLUDEHILL) {
    acceleration.x = force.fx_ex_inner+force.fx_ex_outer;
    acceleration.y = force.fy_ex_inner+force.fy_ex_outer;
    acceleration.z = force.fz_ex_inner+force.fz_ex_outer;
  } 
  else {
    acceleration.x = force.fx_inner+force.fx_outer;
    acceleration.y = force.fy_inner+force.fy_outer;
    acceleration.z = force.fz_inner+force.fz_outer;
  }
  return acceleration;
}

void AdvanceSystemFromDisk(real dt) {
  int NbPlanets, k;
  Point gamma;
  real x, y, z;
  real r, m, smoothing;
  NbPlanets = Sys->nb;
  for (k = 0; k < NbPlanets; k++) {
    if (Sys->FeelDisk[k] == YES) {
      m = Sys->mass[k];
      x = Sys->x[k];
      y = Sys->y[k];
      z = Sys->z[k];
      r = sqrt(x*x + y*y + z*z);
      if (ROCHESMOOTHING != 0)
	smoothing = r*pow(m/3./MSTAR,1./3.)*ROCHESMOOTHING;
      else
	smoothing = ASPECTRATIO*pow(r/R0,FLARINGINDEX)*r*THICKNESSSMOOTHING;
      gamma = ComputeAccel (x, y, z, smoothing, m);
      Sys->vx[k] += dt * gamma.x;
      Sys->vy[k] += dt * gamma.y;
      Sys->vz[k] += dt * gamma.z;
#ifdef GASINDIRECTTERM
      Sys->vx[k] += dt * IndirectTerm.x;
      Sys->vy[k] += dt * IndirectTerm.y;
      Sys->vz[k] += dt * IndirectTerm.z;
#endif
    }
  }
}

OrbitalElements SV2OE (StateVector v, real m) {
  real x,y,z,vx,vy,vz;
  real Ax, Ay, Az, h, h2, inc, e;
  real d, hx, hy, hz, a, E, M, V;
  real hhor, per, an;//Ascending node
  OrbitalElements o;
  x = v.x;
  y = v.y;
  z = v.z;
  vx = v.vx;
  vy = v.vy;
  vz = v.vz;

  d = sqrt(x*x+y*y+z*z);
  
  hx   = y*vz - z*vy;
  hy   = z*vx - x*vz;
  hz   = x*vy - y*vx;
  hhor = sqrt(hx*hx + hy*hy);

  h2  = hx*hx + hy*hy + hz*hz;
  h   = sqrt(h2);
  o.i = inc = asin(hhor/h);

  Ax = vy*hz-vz*hy - G*m*x/d; // v x h - ri/abs(r);
  Ay = vz*hx-vx*hz - G*m*y/d;
  Az = vx*hy-vy*hx - G*m*z/d;

  o.e = e = sqrt(Ax*Ax+Ay*Ay+Az*Az)/(G*m); //Laplace-Runge-Lenz vector
  o.a = a = h*h/(G*m*(1.-e*e));

  //Eccentric anomaly
  if (e != 0.0) {
    E = acos((1.0-d/a)/e); //E evaluated as such is between 0 and PI
  } else {
    E = 0.0;
  }
  if (x*vx+y*vy+z*vz < 0) E= -E; //Planet goes toward central object,
  //hence on its way from aphelion to perihelion (E < 0)

  if (isnan(E)) {
    if (d < a) 
      E = 0.0;
    else
      E = M_PI;
  }

  o.M = M = E-e*sin(E);
  o.E = E;

  //V: true anomaly
  if (e > 1.e-14) {
    V = acos ((a*(1.0-e*e)/d-1.0)/e);
  } else {
    V = 0.0;
  }
  if (E < 0.0) V = -V;

  o.ta = V;
  
  if (fabs(o.i) > 1e-5) {
    an = atan2(hy,hx)+M_PI*.5; //Independently of sign of (hz)
    if (an > 2.0*M_PI) an -= 2.0*M_PI;
  } else {
    an = 0.0;//Line of nodes not determined ==> defaults to x axis
  }

  o.an = an;

  // Argument of periapsis
  per = acos((Ax*cos(an)+Ay*sin(an))/sqrt(Ax*Ax+Ay*Ay+Az*Az));
  if ((-hz*sin(an)*Ax+hz*cos(an)*Ay+(hx*sin(an)-hy*cos(an))*Az) < 0.0)
    per = 2.0*M_PI-per;
  o.per = per;
  if (Ax*Ax+Ay*Ay > 0.0)
    o.Perihelion_Phi = atan2(Ay,Ax);
  else
    o.Perihelion_Phi = atan2(y,x);
  return o;
}

void FindOrbitalElements (v,m,n)
     StateVector v;
     real m;
     int n;
{
  FILE *output;
  char name[256];
  OrbitalElements o;
  if (CPU_Rank) return;
  sprintf (name, "%sorbit%d.dat", OUTPUTDIR, n);
  output = fopen_prs (name, "a");
  o = SV2OE (v,m);
 
  fprintf (output, "%.12g\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g", \
	   PhysicalTime, o.e, o.a, o.M, o.ta, o.per, XAxisRotationAngle);
  fprintf (output, "\t%.12g\t%.12g\t%.12g\n", o.i, o.an, o.Perihelion_Phi);
  fclose (output);
}

void SolveOrbits (sys)
     PlanetarySystem *sys;
{
  int i, n;
  StateVector v;
  n = sys->nb;
  for (i = 0; i < n; i++) {
    v.x = sys->x[i];
    v.y = sys->y[i];
    v.z = sys->z[i];
    v.vx = sys->vx[i];
    v.vy = sys->vy[i];
    v.vz = sys->vz[i];
    FindOrbitalElements (v,MSTAR+sys->mass[i],i);
  }
} 
