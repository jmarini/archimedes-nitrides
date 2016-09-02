/* saveoutput2dgnuplot.h -- This file is part of Archimedes release 0.1.0.
   Archimedes is a simulator for Submicron 2D III-V semiconductor
   Devices. It implements the Monte Carlo method and Hybrid MEP model
   for the simulation of the semiclassical Boltzmann equation for both
   electrons and holes. It also includes the quantum effects by means
   of effective potential method. It is now able to simulate applied
   magnetic fields along with self consistent Faraday equation.

   Copyright (C) 2004, 2005, 2006, 2007 Jean Michel Sellier <sellier@dmi.unict.it>

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 3, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


// ######################################################
// Created on 10 Sep.2004, Siracusa (Italy), J.M.Sellier
// Last modif. : 10 Sep.2007, Siracusa (Italy), J.M.Sellier
// ######################################################

// Here we save all the macroscopic variables
// related to the electrons.
// The format of saving is compatible with GNUPLOT
// and is as follows (constituted of 3 columns):
// x1 y1 u(x1,y1)
// x1 y2 u(x1,y2)
// ...
// xn yn u(xn,yn)
// Actually this file format is completly compatible
// with the "xyz" format of xd3d package. It is higly
// recomended that you use this program in order to
// obtain a very good scientifical visualisation of
// the solutions. In this way the post-processing become
// very easy and clear, so let use this program! It
// is free and you download it from :
// www.gnu.org ---> Free Directory ---> Science ---> Scientific visualisation

// INPUT : je (integer). It is the integer which will be the index of
// the solution. In this way you can easily create a succession of
// solutions like this:
// density000.xyz, density001.xyz, density002.xyz, ..., density010.xyz
// and so on, in order to create movies or animated gif.

void SaveOutput2DGNUPLOT(int je) {
    char s[150];
    FILE *fp;
    FILE *up;
    FILE *vp;
    FILE *lp;
    FILE *lxp;
    FILE *lyp;
    FILE *ep;
    FILE *qp;
    FILE *fM;
    FILE *vo;
    register int i, j;
    int nx = g_mesh->nx,
        ny = g_mesh->ny;
    real dx = g_mesh->dx,
         dy = g_mesh->dy;

    sprintf(s, "density%03d.xyz", je);
    fp = fopen(s, "w");
    sprintf(s, "x_velocity%03d.xyz", je);
    up = fopen(s, "w");
    sprintf(s, "y_velocity%03d.xyz", je);
    vp = fopen(s, "w");
    sprintf(s, "energy%03d.xyz", je);
    ep = fopen(s, "w");
    sprintf(s, "potential%03d.xyz", je);
    lp = fopen(s, "w");
    sprintf(s, "magnetic_field%03d.xyz", je);
    fM = fopen(s,"w");
    sprintf(s, "x_Efield%03d.xyz", je);
    lxp = fopen(s, "w");
    sprintf(s, "y_Efield%03d.xyz", je);
    lyp = fopen(s, "w");
    sprintf(s, "quantum_potential%03d.xyz", je);
    qp = fopen(s, "w");
    sprintf(s, "particle_info%03d.xyz", je);
    vo = fopen(s, "w");

// Save the Monte Carlo results
 if(g_config->simulation_model==MCE || g_config->simulation_model==MCEH){
// Electron Density Output
// =======================
  for(j=1;j<=ny+1;j++){
    for(i=1;i<=nx+1;i++)
      fprintf(fp,"%g %g %g\n",1.e6*(i-1.)*dx,1.e6*(j-1.)*dy,g_mesh->info[i][j].e.density);
    fprintf(fp,"\n");
  }
// X-component of electronic velocity output
// =========================================
  for(j=1;j<=ny+1;j++){
    for(i=1;i<=nx+1;i++)
      fprintf(up,"%g %g %g\n",1.e6*(i-1.)*dx,1.e6*(j-1.)*dy,
              moving_average[i][j][2]);
    fprintf(up,"\n");
  }
// Y-component of electronic velocity output
// =========================================
  for (j=1;j<=ny+1;j++){
    for(i=1;i<=nx+1;i++)
      fprintf(vp,"%g %g %g\n",1.e6*(i-1.)*dx,1.e6*(j-1.)*dy,
              moving_average[i][j][3]);
    fprintf(vp,"\n");
  }
// Electrostatic Potential
// =======================
  for(j=1;j<=ny+1;j++){
    for(i=1;i<=nx+1;i++)
      fprintf(lp,"%g %g %g\n",1.e6*(i-1.)*dx,1.e6*(j-1.)*dy,g_mesh->info[i][j].potential);
    fprintf(lp,"\n");
  }
// Magnetic Field
// ==============
  for(j=1;j<=ny+1;j++){
    for(i=1;i<=nx+1;i++)
      fprintf(fM,"%g %g %g\n",1.e6*(i-1.)*dx,1.e6*(j-1.)*dy,g_mesh->info[i][j].magnetic_field);
    fprintf(fM,"\n");
  }
// X-component of electric field
// =============================
  for(j=1;j<=ny+1;j++){
    for(i=1;i<=nx+1;i++)
      fprintf(lxp,"%g %g %g\n",1.e6*(i-1.)*dx,1.e6*(j-1.)*dy,g_mesh->info[i][j].efield.x);
    fprintf(lxp,"\n");
  }
// Y-component of electric field
// =============================
  for(j=1;j<=ny+1;j++){
    for(i=1;i<=nx+1;i++)
      fprintf(lyp,"%g %g %g\n",1.e6*(i-1.)*dx,1.e6*(j-1.)*dy,g_mesh->info[i][j].efield.y);
    fprintf(lyp,"\n");
  }
// Electron Energy (in eV)
// ===============
  for(j=1;j<=ny+1;j++){
    for(i=1;i<=nx+1;i++)
      fprintf(ep,"%g %g %g\n",1.e6*(i-1.)*dx,1.e6*(j-1.)*dy,
              moving_average[i][j][4]);
    fprintf(ep,"\n");
  }
// Quantum Effective Potential
// ===========================
  for(j=1;j<=ny+1;j++){
    for(i=1;i<=nx+1;i++)
      fprintf(qp,"%g %g %g\n",
          1.e6*(i-1.)*dx,1.e6*(j-1.)*dy,u2d[i][j][0]);
    fprintf(qp,"\n");
  }

  // Particle Information
  // ====================
  fprintf(vo, "index id valley kx ky kz energy t x y i j vx vy\n");
  for(i = 1; i <= g_config->num_particles; ++i) {
      // print all particle info
      particle_info_t *info = &particle_info[i];
      fprintf(vo, "%d %lld %d %g %g %g %g %g %g %g %d %d %g %g\n",
              i,
              info->id,
              info->valley,
              info->kx,
              info->ky,
              info->kz,
              info->energy,
              info->t,
              info->x,
              info->y,
              info->i,
              info->j,
              info->vx,
              info->vy);
  }

// Closure of output files
// =======================
   fclose(fp);
   fclose(up);
   fclose(vp);
   fclose(lp);
   fclose(fM);
   fclose(lxp);
   fclose(lyp);
   fclose(ep);
   fclose(qp);
   fclose(vo);
  }
// Save the Hybrid MEP results
 if(g_config->simulation_model==MEPE || g_config->simulation_model==MEPEH){
// Electron Density Output
// =======================
  for(j=1;j<=ny+1;j++){
    for(i=1;i<=nx+1;i++)
      fprintf(fp,"%g %g %g\n",1.e6*(i-1.)*dx,
              1.e6*(j-1.)*dy,u2d[i+2][j+2][1]);
    fprintf(fp,"\n");
  }
// X-component of electronic velocity output
// =========================================
  for(j=1;j<=ny+1;j++){
    for(i=1;i<=nx+1;i++)
      fprintf(up,"%g %g %g\n",1.e6*(i-1.)*dx,1.e6*(j-1.)*dy,
              u2d[i+2][j+2][2]/u2d[i+2][j+2][1]);
    fprintf(up,"\n");
  }
// Y-component of electronic velocity output
// =========================================
  for (j=1;j<=ny+1;j++){
    for(i=1;i<=nx+1;i++)
      fprintf(vp,"%g %g %g\n",1.e6*(i-1.)*dx,1.e6*(j-1.)*dy,
              u2d[i+2][j+2][3]/u2d[i+2][j+2][1]);
    fprintf(vp,"\n");
  }
// Electrostatic Potential
// =======================
  for(j=1;j<=ny+1;j++){
    for(i=1;i<=nx+1;i++)
      fprintf(lp,"%g %g %g\n",1.e6*(i-1.)*dx,1.e6*(j-1.)*dy,PSI[i][j]);
    fprintf(lp,"\n");
  }
// X-component of electric field
// =============================
  for(j=1;j<=ny+1;j++){
    for(i=1;i<=nx+1;i++)
      fprintf(lxp,"%g %g %g\n",1.e6*(i-1.)*dx,1.e6*(j-1.)*dy,E[i][j][0]);
    fprintf(lxp,"\n");
  }
// Y-component of electric field
// =============================
  for(j=1;j<=ny+1;j++){
    for(i=1;i<=nx+1;i++)
      fprintf(lyp,"%g %g %g\n",1.e6*(i-1.)*dx,1.e6*(j-1.)*dy,E[i][j][1]);
    fprintf(lyp,"\n");
  }
// Electron Energy (in eV)
// ===============
  for(j=1;j<=ny+1;j++){
    for(i=1;i<=nx+1;i++)
      fprintf(ep,"%g %g %g\n",1.e6*(i-1.)*dx,1.e6*(j-1.)*dy,
              u2d[i+2][j+2][4]/u2d[i+2][j+2][1]/Q);
    fprintf(ep,"\n");
  }
// Quantum Effective Potential
// ===========================
  for(j=1;j<=ny+1;j++){
    for(i=1;i<=nx+1;i++)
      fprintf(qp,"%g %g %g\n",
          1.e6*(i-1.)*dx,1.e6*(j-1.)*dy,u2d[i][j][0]);
    fprintf(qp,"\n");
  }

// Closure of output files
// =======================
   fclose(fp);
   fclose(up);
   fclose(vp);
   fclose(lp);
   fclose(lxp);
   fclose(lyp);
   fclose(ep);
   fclose(qp);
  }

}
