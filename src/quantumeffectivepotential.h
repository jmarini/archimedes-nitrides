/* quantumeffectivepotential.h -- This file is part of GNU archimedes.

   Archimedes is a simulator for Submicron and Nanoscaled
   2D III-V Semiconductor Devices.

   Copyright (C) 2004-2011 Jean Michel D. Sellier
   <jeanmichel.sellier@gmail.com>
   <jsellier@purdue.edu>

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
// Created on 28 aug.2004, Siracusa, J.M.Sellier
// Last modif. : 06 Sep.2011, Carry le Rouet, J.M.Sellier
// ######################################################

// For details see:

// QEP_FULL
// ========
// D.K. Ferry, "Effective potentials and the onset of quantization in ultrasmall MOSFETs",
// Superlattices and Microstructures, Vol. 28, No. 5/6, 2000, Equation (2)

// DENSITY_GRADIENT
// ================
// "Schroedinger approach and Density Gradient Model for Quantum effects modeling",
// Simulation Standard, Vol.14, Num.2, February 2004, several SILVACO authors.
// and eventually
// D.K. Ferry, "Effective potentials and the onset of quantization in ultrasmall MOSFETs",
// Superlattices and Microstructures, Vol. 28, No. 5/6, 2000, Equation (10)

// QEP_BOHM
// ========
// D.Bohm, Phys.Rev.85, 166 (1952); 85, 180 (1952).

// QEP_WEIGHTED_BOHM
// =================
// "Effective Bohm Quantum Potential for device simulators based on drift-diffusion and energy transport",
// Simulation os Semiconductor Processes and Devices 2004, vol. 2004, pp.275-278, Munich 2004,
// G.Iannacone, G.Curatola,G.Fiori

void quantum_effective_potential(void)
{
 register int i,j,k;
 real small=1.e-24;
 real n_x[NXM+1][NYM+1];
 real n_y[NXM+1][NYM+1];
 real u[NXM+1][NYM+1];
 real p[5];
 real tmp;
 real Qeff[NXM+1][NYM+1];
 real odx2=1./(dx*dx);
 real ody2=1./(dy*dy);
 real C,D,E;
 real A1=5.;
 real A2=10.;
 real A3;
 real B1=10.;
 real B2=30.;
 real B3;
 real a,b,DELTA;

// WEIGHTED_BOHM and BOHM models
// =============================
 if(QEP_MODEL==QEP_CALIBRATED_BOHM || QEP_MODEL==QEP_BOHM){

  if(QEP_MODEL==QEP_BOHM){
   QEP_GAMMA=1.0;
   QEP_ALPHA=0.5;
  }

  C=-0.5*HBAR*HBAR*QEP_GAMMA/Q;
  D=QEP_ALPHA-1.;
  E=QEP_ALPHA;

  // "normalization" of density
  for(i=1;i<=nx+2;i++)
   for(j=1;j<=ny+2;j++){
    n_x[i][j]=u2d[i][j][1]*small;
    n_y[i][j]=u2d[i][j][1]*small;
   }

  // reconstruction of a smoother density
  // to obtain a more realistic derivative
  // of density w.r.t to position.
  for(j=1;j<=ny+1;j++)
   for(i=1;i<=nx+1-5;i+=5){
    int l;
    for(l=0;l<=4;l++) p[l]=n_x[i+l][j];
    // smooth reconstruction
    A3=-(p[0]+p[1]+p[2]+p[3]+p[4]);
    B3=-(p[1]+2.*p[2]+3.*p[3]+4.*p[4]);
    DELTA=A1*B2-A2*B1;
    if(DELTA!=0.0){
     a=(A2*B3-A3*B2)/DELTA;
     b=(A3*B1-A1*B3)/DELTA;
    } else {
     a=0.2*(p[0]+p[1]+p[2]+p[3]+p[4]);
     b=0.;
    }
    // put it back in the "normalised" density array
    for(l=0;l<=4;l++){
     real val=a+b*l;
     if(val>0.) n_x[i+l][j]=val;
     else n_x[i+l][j]=0.;
    }
   }
  for(i=1;i<=nx+1;i++)
   for(j=1;j<=ny+1-5;j+=5){
    int l;
    for(l=0;l<=4;l++) p[l]=n_y[i][j+l];
    // smooth reconstruction
    A3=-(p[0]+p[1]+p[2]+p[3]+p[4]);
    B3=-(p[1]+2.*p[2]+3.*p[3]+4.*p[4]);
    DELTA=A1*B2-A2*B1;
    if(DELTA!=0.0){
     a=(A2*B3-A3*B2)/DELTA;
     b=(A3*B1-A1*B3)/DELTA;
    } else {
     a=0.2*(p[0]+p[1]+p[2]+p[3]+p[4]);
     b=0.;
    }
    // put it back in the "normalised" density array
    for(l=0;l<=4;l++){
     real val=a+b*l;
     if(val>0.) n_y[i][j+l]=val;
     else n_y[i][j+l]=0.;
    }
   }


   // smooth a bit more to avoid derivative discontinuities
  for(k=0;k<32;k++){ // after many tests, this seems to be a good approximation
   for(i=1;i<=nx+1;i++) for(j=1;j<=ny+1;j++) u[i][j]=n_x[i][j];
   for(i=2;i<=nx;i++)
    for(j=2;j<=ny;j++){
     n_x[i][j]=(u[i][j]+u[i+1][j]+u[i-1][j])/3.;
    }
  }
  for(k=0;k<32;k++){ // after many tests, this seems to be a good approximation
   for(i=1;i<=nx+1;i++) for(j=1;j<=ny+1;j++) u[i][j]=n_y[i][j];
   for(i=2;i<=nx;i++)
    for(j=2;j<=ny;j++){
     n_y[i][j]=(u[i][j]+u[i][j+1]+u[i][j-1])/3.;
    }
  }

  // calculation of effective potential
  // for(i=1;i<=nx+1;i++) for(j=1;j<=ny+1;j++) Qeff[i][j]=0.0;
  for(i=3;i<=nx-1;i++)
    for(j=3;j<=ny-1;j++){
     if(n_x[i][j]>0.0){
      tmp=((MSTAR[i_dom[i+1][j]][1]-MSTAR[i_dom[i][j]][1])/(MSTAR[i_dom[i+1][j]][1]*MSTAR[i_dom[i][j]][1]*M))
                *E*pow(n_x[i][j],D)*(n_x[i+1][j]-n_x[i][j])*odx2
               +((MSTAR[i_dom[i][j+1]][1]-MSTAR[i_dom[i][j]][1])/(MSTAR[i_dom[i][j+1]][1]*MSTAR[i_dom[i][j]][1]*M))
                *E*pow(n_y[i][j],D)*(n_y[i][j+1]-n_y[i][j])*ody2
                +((pow(n_x[i+1][j],QEP_ALPHA)-2.*pow(n_x[i][j],QEP_ALPHA)+pow(n_x[i-1][j],QEP_ALPHA))*odx2
                +(pow(n_y[i][j+1],QEP_ALPHA)-2.*pow(n_y[i][j],QEP_ALPHA)+pow(n_y[i][j-1],QEP_ALPHA))*ody2)
               /(MSTAR[i_dom[i][j]][1]*M);
      Qeff[i][j]=C*tmp/pow(0.5*(n_x[i][j]+n_y[i][j]),QEP_ALPHA);
     } else {
      // no quantum effective potential when no particles are present
      Qeff[i][j]=0.0;
     }
    }
  // Boundary Conditions for effective potential
  for(i=3;i<=nx-1;i++) for(j=1;j<3;j++) Qeff[i][j]=Qeff[i][3];
  for(i=3;i<=nx-1;i++) for(j=ny;j<ny+2;j++) Qeff[i][j]=Qeff[i][ny-1];
  for(j=3;j<=ny-1;j++) for(i=1;i<3;i++) Qeff[i][j]=Qeff[3][j];
  for(j=3;j<=ny-1;j++) for(i=nx;i<nx+2;i++) Qeff[i][j]=Qeff[nx-1][j];
 } // end of weighted_bohm and bohm models

// DENSITY GRADIENT
// ================
 if(QEP_MODEL==QEP_DENSITY_GRADIENT){
  for(i=1;i<=nx+1;i++) for(j=1;j<=ny+1;j++) Qeff[i][j]=0.;

  QEP_ALPHA=0.5;

  C=-HBAR*HBAR*QEP_GAMMA/Q/6.;
  D=QEP_ALPHA-1.;
  E=QEP_ALPHA;

  // "normalization" of density
  for(i=1;i<=nx+2;i++)
   for(j=1;j<=ny+2;j++){
    n_x[i][j]=u2d[i][j][1]*small;
    n_y[i][j]=u2d[i][j][1]*small;
   }

  // reconstruction of a smoother density
  // to obtain a more realistic derivative
  // of density w.r.t to position.
  for(j=1;j<=ny+1;j++)
   for(i=1;i<=nx+1-5;i+=5){
    int l;
    for(l=0;l<=4;l++) p[l]=n_x[i+l][j];
    // smooth reconstruction
    A3=-(p[0]+p[1]+p[2]+p[3]+p[4]);
    B3=-(p[1]+2.*p[2]+3.*p[3]+4.*p[4]);
    DELTA=A1*B2-A2*B1;
    if(DELTA!=0.0){
     a=(A2*B3-A3*B2)/DELTA;
     b=(A3*B1-A1*B3)/DELTA;
    } else {
     a=0.2*(p[0]+p[1]+p[2]+p[3]+p[4]);
     b=0.;
    }
    // put it back in the "normalised" density array
    for(l=0;l<=4;l++){
     real val=a+b*l;
     if(val>0.) n_x[i+l][j]=val;
     else n_x[i+l][j]=0.;
    }
   }
  for(i=1;i<=nx+1;i++)
   for(j=1;j<=ny+1-5;j+=5){
    int l;
    for(l=0;l<=4;l++) p[l]=n_y[i][j+l];
    // smooth reconstruction
    A3=-(p[0]+p[1]+p[2]+p[3]+p[4]);
    B3=-(p[1]+2.*p[2]+3.*p[3]+4.*p[4]);
    DELTA=A1*B2-A2*B1;
    if(DELTA!=0.0){
     a=(A2*B3-A3*B2)/DELTA;
     b=(A3*B1-A1*B3)/DELTA;
    } else {
     a=0.2*(p[0]+p[1]+p[2]+p[3]+p[4]);
     b=0.;
    }
    // put it back in the "normalised" density array
    for(l=0;l<=4;l++){
     real val=a+b*l;
     if(val>0.) n_y[i][j+l]=val;
     else n_y[i][j+l]=0.;
    }
   }

   // smooth a bit more to avoid derivative discontinuities
  for(k=0;k<32;k++){ // after many tests, this seems to be a good approximation
   for(i=1;i<=nx+1;i++) for(j=1;j<=ny+1;j++) u[i][j]=n_x[i][j];
   for(i=2;i<=nx;i++)
    for(j=2;j<=ny;j++){
     n_x[i][j]=(u[i][j]+u[i+1][j]+u[i-1][j])/3.;
    }
  }
  for(k=0;k<32;k++){ // after many tests, this seems to be a good approximation
   for(i=1;i<=nx+1;i++) for(j=1;j<=ny+1;j++) u[i][j]=n_y[i][j];
   for(i=2;i<=nx;i++)
    for(j=2;j<=ny;j++){
     n_y[i][j]=(u[i][j]+u[i][j+1]+u[i][j-1])/3.;
    }
  }

  // calculation of effective potential
  // for(i=1;i<=nx+1;i++) for(j=1;j<=ny+1;j++) Qeff[i][j]=0.0;
  for(i=3;i<=nx-1;i++)
    for(j=3;j<=ny-1;j++){
     if(n_x[i][j]>0.0){
      tmp=((pow(n_x[i+1][j],QEP_ALPHA)-2.*pow(n_x[i][j],QEP_ALPHA)+pow(n_x[i-1][j],QEP_ALPHA))*odx2
         +(pow(n_y[i][j+1],QEP_ALPHA)-2.*pow(n_y[i][j],QEP_ALPHA)+pow(n_y[i][j-1],QEP_ALPHA))*ody2)
         /(MSTAR[i_dom[i][j]][1]*M);
      Qeff[i][j]=C*tmp/pow(0.5*(n_x[i][j]+n_y[i][j]),QEP_ALPHA);
     } else {
      // no quantum effective potential when no particles are present
      Qeff[i][j]=0.0;
     }
    }

  // Boundary Conditions for effective potential
  for(i=3;i<=nx-1;i++) for(j=1;j<3;j++) Qeff[i][j]=Qeff[i][3];
  for(i=3;i<=nx-1;i++) for(j=ny;j<ny+2;j++) Qeff[i][j]=Qeff[i][ny-1];
  for(j=3;j<=ny-1;j++) for(i=1;i<3;i++) Qeff[i][j]=Qeff[3][j];
  for(j=3;j<=ny-1;j++) for(i=nx;i<nx+2;i++) Qeff[i][j]=Qeff[nx-1][j];
 } // end of density gradient model

// FULL effective potential
// ========================
 if(QEP_MODEL==QEP_FULL){
//  int l,m;
//  real tmp;
  real odx2=1./(dx*dx);
  real ody2=1./(dy*dy);
  real alpha2;

  for(i=0;i<=nx+2;i++) for(j=0;j<=ny+2;j++) Qeff[i][j]=0.;

  for(i=2;i<=nx;i++)
   for(j=2;j<=ny;j++){
    // electron wavelenght
    alpha2=HBAR*HBAR/(8.*M*MSTAR[i_dom[i][j]][1]*KB*TL);
    // Mac-Laurin series up to second derivative (leading correction term)
    Qeff[i][j]=alpha2*((PSI[i+1][j]-2.*PSI[i][j]+PSI[i-1][j])*odx2
                      +(PSI[i][j+1]-2.*PSI[i][j]+PSI[i][j-1])*ody2);
   }
 } // end of full effective potential model

 // adds the quantum correction to the electrostatic potential
 for(i=1;i<=nx+1;i++) for(j=1;j<=ny+1;j++) u2d[i][j][0]=PSI[i][j]+Qeff[i][j];

 if(MAXIMINI==1){
    real maxi,mini;
// Compute the maximum and minimum of Quantum Effective Potential
    maxi=Qeff[8][8];
    mini=Qeff[8][8];
    for(i=3;i<=nx-1;i++)
      for(j=3;j<=ny-1;j++){
        if(Qeff[i][j]>=maxi) maxi=Qeff[i][j];
        if(Qeff[i][j]<=mini) mini=Qeff[i][j];
       }
    printf("Max. Eff. Pot. = %g\n",maxi);
    printf("Min. Eff. Pot. = %g\n",mini);
 }
}

// ===============================================
