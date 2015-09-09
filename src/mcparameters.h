/* mcparameters.h -- This file is part of Archimedes release 1.1.0.
   Archimedes is a simulator for Submicron 2D III-V semiconductor
   Devices. It implements the Monte Carlo method
   for the simulation of the semiclassical Boltzmann equation for both
   electrons and holes. It also includes the quantum effects by means 
   of effective potential method. It is now able to simulate applied
   magnetic fields along with self consistent Faraday equation.

   Copyright (C) 2004-2011 Jean Michel Sellier <jeanmichel.sellier@gmail.com>
 
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
// Last modif. : 18 Aug.2011, Carry le Rouet, France, J.M.Sellier
// ######################################################

// Definition of the variables needed for taking into
// account scatterings from acoustic and optical,non-polar
// phonons (which are the most relevant scatterings in common semiconductors like Silicon)

void
MCparameters(int Material)
{
 real wo,no,dos,aco,oge[7],oga[7];
 real dos1,dos2,am1,am2,cl,deq,dij;
 real hwe,hwij,wij,we,ne,nij;
 real poe,poa,ope,opa,eqe,eqa,qmin,qmax;
 real initialenergy,sei,finalenergy,sef;
 real eps,epf,ep,bimp,cimp,qd;
 real ak,qq,wk;
 int ie,i;
 real z2=4.;

 ISEED = 38467.;  //  initial value for random number generator

// These definitions are valid for every material
 BKTQ=KB*TL/Q; // in eV
 QH=Q/HBAR;

// Material with 2 valleys
// #######################
 if(NOVALLEY[Material]==2){
// Effective mass for the GAMMA and L-Valley, respectively
  am1=MSTAR[Material][1]*M;
  am2=MSTAR[Material][2]*M;
  eps=EPSR[Material]*EPS0;
  epf=EPF[Material]*EPS0;
  ep=1./(1./epf-1./eps);

// Parameters for Phonon Scattering
  cl=RHO[Material]*pow(UL[Material],2.);
  deq=dij=DTK[Material][0]*Q;
  hwe=hwij=HWO[Material][0];

  SMH[Material][1]=sqrt(2.*am1*Q)/HBAR;
  SMH[Material][2]=sqrt(2.*am2*Q)/HBAR;
  HHM[Material][1]=HBAR*HBAR/(2.*am1*Q);
  HHM[Material][2]=HBAR*HBAR/(2.*am2*Q);
  HM[Material][1]=HBAR/am1;
  HM[Material][2]=HBAR/am2;

  wo=HWO[Material][0]*Q/HBAR;
  wij=hwij*Q/HBAR;
  we=hwe*Q/HBAR;

  no=1./(exp(HWO[Material][0]/BKTQ)-1.);
  nij=1./(exp(hwij/BKTQ)-1.);
  ne=1./(exp(hwe/BKTQ)-1.);

  dos1=pow(sqrt(2.*am1*Q)/HBAR,3.)/pow(2.*PI,2.);
  dos2=pow(sqrt(2.*am2*Q)/HBAR,3.)/pow(2.*PI,2.);

  poe=Q/8./PI/ep*Q*wo*(no+1.);
  poa=poe*no/(1.+no);
  aco=2.*PI*DA[Material]/Q*DA[Material]*BKTQ/HBAR*Q/cl;
  ope=PI*dij/wij*dij/RHO[Material]/Q*(nij+1.);
  opa=ope*nij/(1.+nij);
  eqe=PI*deq/we*deq/RHO[Material]/Q*(ne+1.);
  eqa=eqe*ne/(1.+ne);

// Parameters for impurity scatterings
  cimp=CIMP; // <--- impurity concentration
  qd=sqrt(Q*cimp/BKTQ/eps);
  QD2=qd*qd;
  bimp=2.*PI*cimp*Q*Q/HBAR*Q/eps/eps;

// Calculation of scattering rates
  for(ie=1;ie<=DIME;ie++){
   initialenergy=DE*((real)(ie));
   sei=sqrt(initialenergy);
// GAMMA-valley
// ============
   if(OPTICALPHONONS==ON){
// Polar optical phonon
// Emission - Parabolic
    finalenergy=initialenergy-HWO[Material][0];
    if(finalenergy>0.){
     sef=sqrt(finalenergy);
     qmax=sef+sei;
     qmin=sei-sef;
     SWK[Material][1][1][ie]=poe*SMH[Material][1]*sei/initialenergy/Q*log(qmax/qmin);
    }
    else SWK[Material][1][1][ie]=0.;
// Absorption - Parabolic
    finalenergy=initialenergy+HWO[Material][0];
    sef=sqrt(finalenergy);
    qmax=sef+sei;
    qmin=sef-sei;
    SWK[Material][1][2][ie]=SWK[Material][1][1][ie]
                 +poa*SMH[Material][1]*sei/initialenergy/Q*log(qmax/qmin);
// Emission
    finalenergy=initialenergy-hwij+EMIN[Material][1]-EMIN[Material][2];
    if(finalenergy>0.){
     sef=sqrt(finalenergy*(1.+alphaK[Material][2]*finalenergy));
     SWK[Material][1][3][ie]=SWK[Material][1][2][ie]
                            +z2*ope*sef*dos2*(1.+2.*alphaK[Material][2]*finalenergy);
    }
    else SWK[Material][1][3][ie]=SWK[Material][1][2][ie];
// Absorption
    finalenergy=initialenergy+hwij+EMIN[Material][1]-EMIN[Material][2];
    if(finalenergy>0.){
     sef=sqrt(finalenergy*(1.+alphaK[Material][2]*finalenergy));
     SWK[Material][1][4][ie]=SWK[Material][1][3][ie]
                            +z2*opa*sef*dos2*(1.+2.*alphaK[Material][2]*finalenergy);
    }
    else SWK[Material][1][4][ie]=SWK[Material][1][3][ie];
   }
   else{
    // NO OPTICAL SCATTERING
    SWK[Material][1][1][ie]=0.0;
    SWK[Material][1][2][ie]=0.0;
    SWK[Material][1][3][ie]=0.0;
    SWK[Material][1][4][ie]=0.0;
   }
   if(ACOUSTICPHONONS==ON){
// Acoustic Phonon
    finalenergy=initialenergy;
    sef=sqrt(finalenergy*(1.+alphaK[Material][1]*finalenergy));
    SWK[Material][1][5][ie]=SWK[Material][1][4][ie]
                           +aco*sef*dos1*(1.+2.*alphaK[Material][1]*finalenergy);
   }
   else {
    // NO ACOUSTIC PHONON
    SWK[Material][1][5][ie]=SWK[Material][1][4][ie]+0.0;
   }
// Impurity scattering
   if(IMPURITYPHONONS==ON){
    finalenergy=initialenergy;
    sef=sqrt(finalenergy*(1.+alphaK[Material][1]*finalenergy));
    ak=SMH[Material][1]*sef;
    qq=QD2*(4.*ak*ak+QD2);
    wk=bimp/qq*sef*dos1*(1.+2.*alphaK[Material][1]*finalenergy);
    if(wk>1.e14) wk=1.e14;
    SWK[Material][1][6][ie]=SWK[Material][1][5][ie]+wk;
   }
   else {
    // NO IMPURITY SCATTERING
    SWK[Material][1][6][ie]=SWK[Material][1][5][ie]+0.0;
   }
// L-valley
// ========
   if(OPTICALPHONONS==ON){
// Polar optical phonon
// Emission
    finalenergy=initialenergy-HWO[Material][0];
    if(finalenergy>0.){
      sef=sqrt(finalenergy);
      qmax=sef+sei;
      qmin=sei-sef;
      SWK[Material][2][1][ie]=poe*SMH[Material][2]*sei/initialenergy/Q*log(qmax/qmin);
    }
    else SWK[Material][2][1][ie]=0.;
// Absorption
    finalenergy=initialenergy+HWO[Material][0];
    sef=sqrt(finalenergy);
    qmax=sef+sei;
    qmin=sef-sei;
    SWK[Material][2][2][ie]=SWK[Material][2][1][ie]
                 +poa*SMH[Material][2]*sei/initialenergy/Q*log(qmax/qmin);
// Non-polar optical phonon
// Emission
    finalenergy=initialenergy-hwe;
    if(finalenergy>0.){
      sef=sqrt(finalenergy*(1.+alphaK[Material][2]*finalenergy));
      SWK[Material][2][3][ie]=SWK[Material][2][2][ie]
                   +(z2-1.)*eqe*sef*dos2*(1.+2.*alphaK[Material][2]*finalenergy);
    }
    else SWK[Material][2][3][ie]=SWK[Material][2][2][ie];
// Absorption
    finalenergy=initialenergy+hwe;
    sef=sqrt(finalenergy*(1.+alphaK[Material][2]*finalenergy));
    SWK[Material][2][4][ie]=SWK[Material][2][3][ie]
                 +(z2-1.)*eqa*sef*dos2*(1.+2.*alphaK[Material][2]*finalenergy);

// Emission
    finalenergy=initialenergy-hwij+EMIN[Material][2]-EMIN[Material][1];
    if(finalenergy>0.){
     sef=sqrt(finalenergy*(1.+alphaK[Material][2]*finalenergy));
     SWK[Material][2][5][ie]=SWK[Material][2][4][ie]
                            +ope*sef*dos1*(1.+2.*alphaK[Material][2]*finalenergy);
    }
    else SWK[Material][2][5][ie]=SWK[Material][2][4][ie];
// Absorption
    finalenergy=initialenergy+hwij+EMIN[Material][2]-EMIN[Material][1];
    if(finalenergy>0.){
     sef=sqrt(finalenergy*(1.+alphaK[Material][2]*finalenergy));
     SWK[Material][2][6][ie]=SWK[Material][2][5][ie]
                            +opa*sef*dos1*(1.+2.*alphaK[Material][2]*finalenergy);
    }
    else SWK[Material][2][6][ie]=SWK[Material][2][5][ie];
   }
   else {
      // NO OPTICAL SCATTERING
      SWK[Material][2][1][ie]=0.0;
      SWK[Material][2][2][ie]=0.0;
      SWK[Material][2][3][ie]=0.0;
      SWK[Material][2][4][ie]=0.0;
      SWK[Material][2][5][ie]=0.0;
      SWK[Material][2][6][ie]=0.0;
   }
   if(ACOUSTICPHONONS==ON){
// Acoustic phonon
    finalenergy=initialenergy;
    sef=sqrt(finalenergy*(1.+alphaK[Material][2]*finalenergy));
    SWK[Material][2][7][ie]=SWK[Material][2][6][ie]
                           +aco*sef*dos2*(1.+2.*alphaK[Material][2]*finalenergy);
   }
   else {
      // NO ACOUSTIC SCATTERING
      SWK[Material][2][7][ie]=SWK[Material][2][6][ie]+0.0;
   }
// ==================
// --->  } <--- This is a bug! Removed by J.M.Sellier on 24 dec.2006 (SR)
// (Thanks to Kun-Yuan Xu who helped me to find it!)
// ==================
   if(IMPURITYPHONONS==ON){
// Impurity scattering
    finalenergy=initialenergy;
    sef=sqrt(finalenergy*(1.+alphaK[Material][2]*finalenergy));
    ak=SMH[Material][2]*sef;
    qq=QD2*(4.*ak*ak+QD2);
    wk=bimp/qq*sef*dos2*(1.+2.*alphaK[Material][2]*finalenergy);
    if(wk>1.e14) wk=1.e14;
    SWK[Material][2][8][ie]=SWK[Material][2][7][ie]+wk;
   }
   else {
    // NO IMPURITY SCATTERING
    SWK[Material][2][8][ie]=SWK[Material][2][7][ie]+0.0;
   }
  } // <--- this is the correct place for the symbol } (J.M.Sellier,24/12/2006)
// Evalutation of gamma
  GM[Material]=SWK[Material][1][6][1];
  for(ie=1;ie<=DIME;ie++){
    if(SWK[Material][1][6][ie]>GM[Material]) GM[Material]=SWK[Material][1][6][ie];
    if(SWK[Material][2][8][ie]>GM[Material]) GM[Material]=SWK[Material][2][8][ie];
  }
  printf("GAMMA[%d] = %g \n",Material, GM[Material]);
  for(i=1;i<=6;i++)
    for(ie=1;ie<=DIME;ie++)
      SWK[Material][1][i][ie]/=GM[Material];
  for(i=1;i<=8;i++)
    for(ie=1;ie<=DIME;ie++)
      SWK[Material][2][i][ie]/=GM[Material];
 }
// End of two-valleys materials
// ############################

// Material = one valley
// #####################
 if(NOVALLEY[Material]==1){
  SMH[Material][0]=sqrt(2.*MSTAR[Material][1]*M*Q)/HBAR;
  HHM[Material][0]=HBAR*HBAR/(2.*MSTAR[Material][1]*M*Q);
  HM[Material][0]=HBAR/(MSTAR[Material][1]*M);
// Density of states
  dos=pow((sqrt(2.*MSTAR[Material][1]*M)*sqrt(Q)/HBAR),3.)/(4.*PI*PI);
// constant for the acoustic phonon
  aco=2.*PI*(DA[Material]/Q)*DA[Material]*(BKTQ/HBAR)
     *(Q/(RHO[Material]*UL[Material]*UL[Material]));
// Constants for the 6 Silicon-like non-polar optical phonons
   for(i=1;i<=6;i++){
// i-th Optical Phonon
    oge[i]=0.;
    oga[i]=0.;
    if(ZF[Material][i-1]!=0.){
      wo=HWO[Material][i-1]*Q/HBAR; // frequency of phonon
      no=1./(exp(HWO[Material][i-1]/BKTQ) - 1.); // population of phonons
      oge[i]=ZF[Material][i-1]*PI*(DTK[Material][i-1]*Q/wo)
            *((DTK[Material][i-1]*Q/RHO[Material])/Q)*(no+1.);
      oga[i]=oge[i]*no/(1.+no);
    }
   }
// Calculation of scattering rates
   for(ie=1; ie<=DIME; ++ie) SWK[Material][0][0][ie]=0.;
   for(ie=1; ie<=DIME; ++ie){
    initialenergy=DE*((real) ie);
    sei=sqrt(initialenergy);
    if(OPTICALPHONONS==ON){
// non polar optical phonons
     for(i=1;i<=6;i++){
      finalenergy=initialenergy-HWO[Material][i-1];
      SWK[Material][0][i*2-1][ie]=SWK[Material][0][i*2-2][ie];
      if(finalenergy>0.){
       sef=sqrt(finalenergy*(1.+alphaK[Material][1]*finalenergy));
       SWK[Material][0][i*2-1][ie]=SWK[Material][0][i*2-2][ie]
                     +oge[i]*sef*dos*(1.+2.*alphaK[Material][1]*finalenergy);
      }
      finalenergy=initialenergy+HWO[Material][i-1];
      SWK[Material][0][i*2][ie]=SWK[Material][0][i*2-1][ie];
      if(finalenergy>0.){
       sef=sqrt(finalenergy*(1.+alphaK[Material][1]*finalenergy));
       SWK[Material][0][i*2][ie]=SWK[Material][0][i*2-1][ie]
                   +oga[i]*sef*dos*(1.+2.*alphaK[Material][1]*finalenergy);
      }
     }
    }
    else {
     // NO OPTICAL SCATTERING
     SWK[Material][0][1][ie]=0.0;
     SWK[Material][0][2][ie]=0.0;
     SWK[Material][0][3][ie]=0.0;
     SWK[Material][0][4][ie]=0.0;
     SWK[Material][0][5][ie]=0.0;
     SWK[Material][0][6][ie]=0.0;
     SWK[Material][0][7][ie]=0.0;
     SWK[Material][0][8][ie]=0.0;
     SWK[Material][0][9][ie]=0.0;
     SWK[Material][0][10][ie]=0.0;
     SWK[Material][0][11][ie]=0.0;
     SWK[Material][0][12][ie]=0.0;
    }
    if(ACOUSTICPHONONS==ON){
// Acoustic phonon
    finalenergy=initialenergy;
    sef=sqrt(finalenergy*(1.+alphaK[Material][1]*finalenergy));
    SWK[Material][0][13][ie]=SWK[Material][0][12][ie]
                  +aco*sef*dos*(1.+2.*alphaK[Material][1]*finalenergy);
    }
    else {
    // NO ACOUSTIC PHONONS 
    SWK[Material][0][13][ie]=SWK[Material][0][12][ie]+0.0;
    }
   }
// Evaluation of gamma
  GM[Material]=SWK[Material][0][13][1];
  for(ie=1;ie<=DIME;++ie)
    if(SWK[Material][0][13][ie]>GM[Material]) GM[Material]=SWK[Material][0][13][ie];
  printf("GAMMA[%d] = %g \n",Material, GM[Material]);
  for(ie=1;ie<=DIME;ie++)
    for(i=1;i<=13;i++)
      SWK[Material][0][i][ie]/=GM[Material];
 }
// End of one-valley material
}

// =======================================================
