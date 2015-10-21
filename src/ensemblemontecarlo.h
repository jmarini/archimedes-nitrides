/* ensemblemontecarlo.h -- This file is part of Archimedes release 1.2.0.
   Archimedes is a simulator for Submicron 2D III-V semiconductor
   Devices. It implements the Monte Carlo method
   for the simulation of the semiclassical Boltzmann equation for both
   electrons and holes. It includes some quantum effects by means 
   of effective potential method. It is also able to simulate applied
   magnetic fields along with self consistent Faraday equation.

   Copyright (C) 2004-2011 Jean Michel Sellier
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
// Created on 05 sep.2004, Siracusa, J.M.Sellier
// Last modif. : 25 Aug.2011, Carry le Rouet, France, J.M.Sellier
// ######################################################

// Ensemble Monte Carlo method

void
EMC(void)
{
 long int n=1;
 int i,ni,j,npt[NXM+NYM+1][4];
 real tdt,ti,tau;

 memset(&npt,0,sizeof(npt));
 tdt=TEMPO+DT;

 do{
    // information about particle n is set up in easy access variables
    IV = P[n].valley;
    KX = P[n].kx;
    KY = P[n].ky;
    KZ = P[n].kz;
    TS = P[n].t;
    X  = P[n].x;
    Y  = P[n].y;
    ti = TEMPO;

    // while the particle's time is less than the time for the step...
    while(TS<=tdt){
       tau=TS-ti;          // the dt for the current step
       drift(tau);         // drift for dt
       i=(int)(X/dx)+1;    // update cell position after drift -- i
       j=(int)(Y/dy)+1;    //                                  -- j
       if(i<=1)  i=1;      // make sure variables are in bound -- i
       if(j<=1)  j=1;      //                                  -- j
       if(i>=nx) i=nx;     //                                  -- nx
       if(j>=ny) j=ny;     //                                  -- ny
       scat(i_dom[i][j]);  // scatter particle
       ti=TS;              // update the time
       i=(int)(X/dx)+1;    // update cell position after scatter -- i
       j=(int)(Y/dy)+1;    //                                    -- j
       if(i<=1)  i=1;      // make sure variables are in bound -- i
       if(j<=1)  j=1;      //                                  -- j
       if(i>=nx) i=nx;     //                                  -- nx
       if(j>=ny) j=ny;     //                                  -- ny
       TS=ti-log(rnd())/GM[i_dom[i][j]];
    }
    tau=tdt-ti;            // calculate unused time in step
    drift(tau);            // drift for unused time in step

    // check if a particle is going out from the right edge of the device
    if(IV!=9){
     i=(int)(X/dx+1.5);
     j=(int)(Y/dy+1.5);
      if(i>=nx+1 && (EDGE[1][j][0]==1 || EDGE[1][j][0]==2)){
        IV=9;
        if(npt[j][1]<(NP1/2) && j>1 && j<ny+1){
          npt[j][1]++;
          IV=1;
        }
        else if(npt[j][1]<(NP1/4) && (j<=1 || j>=ny+1)){
          npt[j][1]++;
          IV=1;
        }
      }
    }
    // check if a particle is going out from the left edge of the device
    if(IV!=9){
     i=(int)(X/dx+1.5);
     j=(int)(Y/dy+1.5);
      if(i<=1 && (EDGE[3][j][0]==1 || EDGE[3][j][0]==2)){
        IV=9;
        if(npt[j][3]<(NP1/2) && j>1 && j<ny+1){
           npt[j][3]++;
           IV=1;
        }
        else if(npt[j][3]<(NP1/4) && (j<=1 || j>=ny+1)){
          npt[j][3]++;
          IV=1;
        }
      }
    }
    // check if a particle is going out from the bottom edge of the device
    if(IV!=9){
     i=(int)(X/dx+1.5);
     j=(int)(Y/dy+1.5);
      if(j<=1 && (EDGE[0][i][0]==1 || EDGE[0][i][0]==2)){
        IV=9;
        if(npt[i][0]<(NP1/2) && (i>1 || i<nx+1)){
           npt[i][0]++;
           IV=1;
        }
        if(npt[i][0]<(NP1/4) && (i<=1 || i>=nx+1)){
           npt[i][0]++;
           IV=1;
        }
      }
    }
    // check if a particle is going out from the upper edge of the device
    if(IV!=9){
     i=(int)(X/dx+1.5);
     j=(int)(Y/dy+1.5);
      if(j>=ny+1 && (EDGE[2][i][0]==1 || EDGE[2][i][0]==2)){
        IV=9;
        if(npt[i][2]<(NP1/2) && (i>1 || i<nx+1)){
           npt[i][2]++;
           IV=1;
        }
        if(npt[i][2]<(NP1/4) && (i<=1 || i>=nx+1)){
           npt[i][2]++;
           IV=1;
        }
      }
    }


// ==============================================
     if(IV!=9){
       P[n].valley = IV;
       P[n].kx = KX;
       P[n].ky = KY;
       P[n].kz = KZ;
       P[n].t  = TS;
       P[n].x  = X;
       P[n].y  = Y;
       n++;
     }
// if IV=9 then the super-particle has been eliminated
     if(IV==9){
       P[n] = P[INUM];
       INUM--;
     }
  }while(n<INUM);

// create particles at ohmic contacts of the bottom edge
  for(i=1;i<=nx+1;i++){
//printf("\n\n\n********************* %d ************************\n\n\n",npt[i][0]);
   if(EDGE[0][i][0]==2){
    ni=(NP1/2)-npt[i][0];
    if(i==1 || i==nx+1) ni=NP1/4-npt[i][0];
    if(ni>0){
      for(j=1;j<=ni;j++){
        n=INUM+j;
        creation(i,TEMPO,0);
        P[n].valley = IV;
        P[n].kx = KX;
        P[n].ky = KY;
        P[n].kz = KZ;
        P[n].t  = TS;
        P[n].x  = X;
        P[n].y  = Y;
      }
      INUM += ni;
    }
   }
  }

// create particles at ohmic contacts of the upper edge
  for(i=1;i<=nx+1;i++){
   if(EDGE[2][i][0]==2){
    ni=(NP1/2)-npt[i][2];
    if(i==1 || i==nx+1) ni=NP1/4-npt[i][2];
    if(ni>0){
      for(j=1;j<=ni;j++){
        n=INUM+j;
        creation(i,TEMPO,2);
        P[n].valley = IV;
        P[n].kx = KX;
        P[n].ky = KY;
        P[n].kz = KZ;
        P[n].t  = TS;
        P[n].x  = X;
        P[n].y  = Y;
      }
      INUM += ni;
    }
   }
  }
// create particles at ohmic contacts of the right edge
  for(i=1;i<=ny+1;i++){
   if(EDGE[1][i][0]==2){
    ni=(NP1/2)-npt[i][1];
    if(i==1 || i==ny+1) ni=NP1/4-npt[i][1];
    if(ni>0){
      for(j=1;j<=ni;j++){
        n=INUM+j;
        creation(i,TEMPO,1);
        P[n].valley = IV;
        P[n].kx = KX;
        P[n].ky = KY;
        P[n].kz = KZ;
        P[n].t  = TS;
        P[n].x  = X;
        P[n].y  = Y;
      }
      INUM += ni;
    }
   }
  }
// create particles at ohmic contacts of the left edge
  for(i=1;i<=ny+1;i++){
   if(EDGE[3][i][0]==2){
    ni=(NP1/2)-npt[i][3];
    if(i==1 || i==ny+1) ni=NP1/4-npt[i][3];
    if(ni>0){
      for(j=1;j<=ni;j++){
        n=INUM+j;
        creation(i,TEMPO,3);
        P[n].valley = IV;
        P[n].kx = KX;
        P[n].ky = KY;
        P[n].kz = KZ;
        P[n].t  = TS;
        P[n].x  = X;
        P[n].y  = Y;
      }
      INUM += ni;
    }
   }
  }

 printf("Actual number of electron super-particles = %d\n",INUM);
 if(INUM>NPMAX){
   printf("%s: too big actual number of particles\n",progname);
   exit(EXIT_FAILURE);
 }
}

// ============================================================
