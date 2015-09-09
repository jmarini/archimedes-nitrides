/* scattering.h -- This file is part of Archimedes release 1.2.0
   Archimedes is a simulator for Submicron 2D III-V semiconductor
   Devices. It implements the Monte Carlo method
   for the simulation of the semiclassical Boltzmann equation for both
   electrons and holes. It also includes some quantum effects by means 
   of the effective potential method. It is also able to simulate applied
   magnetic fields along with self consistent Faraday equation.

   Copyright (C) 2004-2011 Jean Michel Sellier
   <jeanmichel.sellier@gmail.com>
 
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
// Created on 06 sep.2004, Siracusa, J.M.Sellier
// Last modif. : 31 Aug.2011, Carry le Rouet, J.M.Sellier
// ######################################################

// calculation of the various scattering process
// From version 0.0.7 on, we take into account
// acoustic phonons scattering and all the non-polar optical phonons
// scattering which are the most relevant in one-valley materials,
// and two-valleys materials.
// From version 1.1.0 on, the scattering effects can be excluded
// to simulate ballistic transport.

void
scat(int Material)
{
 int j=0;
 register int i,ie=0;
 real ksquared,thesquareroot,superparticle_energy;
 real r1,finalenergy=0.,finalk,cosinus,sinus,fai;
 real f,cb,cf,sf,skk,a11,a12,a13,a21,a22,a23,a32,a33,x1,x2,x3,sb,r2;
 real ki,kf,cs,sn;

  if(IV==9) return;

// One-valley material Scattering Selection
// ########################################
 if(NOVALLEY[Material]==1){
  ksquared=KX*KX+KY*KY+KZ*KZ;
  if(CONDUCTION_BAND==FULL){
   real k,k2,k4;
   k=sqrt(ksquared)*0.5/PI*1.e-12;
   k2=k*k;
   k4=k2*k2;
   k=sqrt(k2);
   // periodicity on reciprocal lattice
   superparticle_energy=CB_FULL[Material][0]*k4*k4*k2
                       +CB_FULL[Material][1]*k4*k4*k
                       +CB_FULL[Material][2]*k4*k4
                       +CB_FULL[Material][3]*k4*k2*k
                       +CB_FULL[Material][4]*k4*k2
                       +CB_FULL[Material][5]*k4*k
                       +CB_FULL[Material][6]*k4
                       +CB_FULL[Material][7]*k2*k
                       +CB_FULL[Material][8]*k2
                       +CB_FULL[Material][9]*k
                       +CB_FULL[Material][10]; // in eV
  }
  if(CONDUCTION_BAND==KANE){
   thesquareroot=sqrt(1.+4.*alphaK[Material][1]*HHM[Material][0]*ksquared);
   superparticle_energy=(thesquareroot-1.)/(2.*alphaK[Material][1]);
  }
  if(CONDUCTION_BAND==PARABOLIC){
   superparticle_energy=HHM[Material][0]*ksquared; // in eV
  }
  if(superparticle_energy<=0.) return;
  ie=((int)(superparticle_energy/DE))+1;
  if(ie>DIME) ie=DIME;

// Selection of scattering process
  r1 = rnd();
// Non-Polar optical phonons
  for(i=1;i<=6;i++){
// Emission of an optical phonon
    if(r1<=SWK[Material][0][i*2-1][ie] && j==0){
       finalenergy=superparticle_energy-HWO[Material][i-1];
       if(finalenergy<=0.) return;
       j=1;
    }
// Absorbation of an optical phonon
    if((r1<=SWK[Material][0][i*2][ie]) && j==0){
       finalenergy=superparticle_energy+HWO[Material][i-1];
       if(finalenergy<=0.) return;
       j=1;
    }
  }
// Acoustic phonon
  if((r1<=SWK[Material][0][13][ie]) && j==0){
     finalenergy=superparticle_energy;
     if(finalenergy<=0.) return;
     j=1;
  }
  if((finalenergy<=0.) || j==0) return;

// determination of the final states
  if(CONDUCTION_BAND==FULL){
   // look for a final k
   // bisection algorithm, probably not best but at least something to start from..
    real x=0.0,a=0.,b=1./LATTCONST[Material]*1.e-12,a2;
    // real del=0.001
    real a4,x4;
    real fx,fa;
    int k=0;
    for(k=0;k<=512;k++){
     x=0.5*(a+b);
     x2=x*x;
     x4=x2*x2;
     a2=a*a;
     a4=a2*a2;
     fa=(CB_FULL[Material][0]*a4*a4*a2
        +CB_FULL[Material][1]*a4*a4*a
        +CB_FULL[Material][2]*a4*a4
        +CB_FULL[Material][3]*a4*a2*a
        +CB_FULL[Material][4]*a4*a2
        +CB_FULL[Material][5]*a4*a
        +CB_FULL[Material][6]*a4
        +CB_FULL[Material][7]*a2*a
        +CB_FULL[Material][8]*a2
        +CB_FULL[Material][9]*a
        +CB_FULL[Material][10])-finalenergy;
     fx=(CB_FULL[Material][0]*x4*x4*x2
        +CB_FULL[Material][1]*x4*x4*x
        +CB_FULL[Material][2]*x4*x4
        +CB_FULL[Material][3]*x4*x2*x
        +CB_FULL[Material][4]*x4*x2
        +CB_FULL[Material][5]*x4*x
        +CB_FULL[Material][6]*x4
        +CB_FULL[Material][7]*x2*x
        +CB_FULL[Material][8]*x2
        +CB_FULL[Material][9]*x
        +CB_FULL[Material][10])-finalenergy;
     if((fa*fx)<0.) b=x;
     else a=x;
    }
   finalk=x*1.e12*2.*PI;
  }
  if(CONDUCTION_BAND==KANE) finalk = SMH[Material][0]*sqrt(finalenergy*(1.+alphaK[Material][1]*finalenergy));
  if(CONDUCTION_BAND==PARABOLIC) finalk=SMH[Material][0]*sqrt(finalenergy);
  cosinus = 1.-2.*rnd();
  sinus = sqrt(1.-cosinus*cosinus);
  fai = 2.*PI*rnd();
  KX = finalk*cosinus;
  KY = finalk*sinus*cos(fai);
  KZ = finalk*sinus*sin(fai);
  return;
 }

// Two-valleys material Scattering Selection
// #########################################
 if(NOVALLEY[Material]==2){
   ksquared=KX*KX+KY*KY+KZ*KZ;
   ki=sqrt(ksquared);
   if(CONDUCTION_BAND==KANE){
    thesquareroot=sqrt(1.+4.*alphaK[Material][IV]*HHM[Material][IV]*ksquared);
    superparticle_energy=(thesquareroot-1.)/(2.*alphaK[Material][IV]);
   }
   if(CONDUCTION_BAND==PARABOLIC){
    superparticle_energy=HHM[Material][IV]*ksquared;
   }
   if(superparticle_energy<=0.) return;
   ie=((int)(superparticle_energy/DE))+1;
   if(ie>DIME) ie=DIME;

// Selection of scattering process in the GAMMA-Valley
// ===================================================
  if(IV==1){
    r1 = rnd();
// Non-Polar optical phonons
// Emission of an optical phonon
      if(r1<=SWK[Material][1][1][ie] && j==0){
         finalenergy=superparticle_energy-HWO[Material][0];
         if(finalenergy<=0.) return;
//         if(finalenergy>DIME*DE) finalenergy=BKTQ;
         j=1;
// linea 20
       if(CONDUCTION_BAND==KANE) kf=SMH[Material][IV]*sqrt(finalenergy*(1.+alphaK[Material][IV]*finalenergy));
       if(CONDUCTION_BAND==PARABOLIC) kf=SMH[Material][IV]*sqrt(finalenergy);
       f=2.*ki*kf/(ki-kf)/(ki-kf);
       if(f<=0.) return;
       cb=(1.+f-pow(1.+2.*f,rnd()))/f;
// linea 30 -- determination of the final states
       sb=sqrt(1.-cb*cb);
       fai=2.*PI*rnd();
       cf=cos(fai);
       sf=sin(fai);
       skk=sqrt(KX*KX+KY*KY);
       a11=KY/skk;
       a12=KX*KZ/skk/ki;
       a13=KX/ki;
       a21=-KX/skk;
       a22=KY*KZ/skk/ki;
       a23=KY/ki;
       a32=-skk/ki;
       a33=KZ/ki;
       x1=kf*sb*cf;
       x2=kf*sb*sf;
       x3=kf*cb;
       KX=a11*x1+a12*x2+a13*x3;
       KY=a21*x1+a22*x2+a23*x3;
       KZ=a32*x2+a33*x3;
       return;
      }
// Absorption of an optical phonon
      if((r1<=SWK[Material][1][2][ie]) && j==0){
         finalenergy=superparticle_energy+HWO[Material][0];
//         if(finalenergy<0.) return;
//         if(finalenergy>DIME*DE) finalenergy=BKTQ;
         j=1;
// linea 20
       if(CONDUCTION_BAND==KANE) kf=SMH[Material][IV]*sqrt(finalenergy*(1.+alphaK[Material][IV]*finalenergy));
       if(CONDUCTION_BAND==PARABOLIC) kf=SMH[Material][IV]*sqrt(finalenergy);
       f=2.*ki*kf/(ki-kf)/(ki-kf);
       if(f<=0.) return;
       cb=(1.+f-pow((1.+2.*f),rnd()))/f;
// linea 30 -- determination of the final states
       sb=sqrt(1.-cb*cb);
       fai=2.*PI*rnd();
       cf=cos(fai);
       sf=sin(fai);
       skk=sqrt(KX*KX+KY*KY);
       a11=KY/skk;
       a12=KX*KZ/skk/ki;
       a13=KX/ki;
       a21=-KX/skk;
       a22=KY*KZ/skk/ki;
       a23=KY/ki;
       a32=-skk/ki;
       a33=KZ/ki;
       x1=kf*sb*cf;
       x2=kf*sb*sf;
       x3=kf*cb;
       KX=a11*x1+a12*x2+a13*x3;
       KY=a21*x1+a22*x2+a23*x3;
       KZ=a32*x2+a33*x3;
       return;
      }
// Emission
      if((r1<=SWK[Material][1][3][ie]) && j==0){
         finalenergy=superparticle_energy-HWO[Material][0]+EMIN[Material][1]-EMIN[Material][2];
         if(finalenergy<=0.) return;
//         if(finalenergy>DIME*DE) finalenergy=BKTQ;
         IV=2;
         j=1;
// determination of the final states
       if(CONDUCTION_BAND==KANE) kf = SMH[Material][IV]*sqrt(finalenergy*(1.+alphaK[Material][IV]*finalenergy));
       if(CONDUCTION_BAND==PARABOLIC) kf=SMH[Material][IV]*sqrt(finalenergy);
       cs = 1.-2.*rnd();
       sn = sqrt(1.-cs*cs);
       fai = 2.*PI*rnd();
       KX = kf*cs;
       KY = kf*sn*cos(fai);
       KZ = kf*sn*sin(fai);
       return;
      }
// Absorption
      if((r1<=SWK[Material][1][4][ie]) && j==0){
         finalenergy=superparticle_energy+HWO[Material][0]+EMIN[Material][1]-EMIN[Material][2];
         if(finalenergy<=0.) return;
//         if(finalenergy>DIME*DE) finalenergy=BKTQ;
         IV=2;
         j=1;
// determination of the final states
       if(CONDUCTION_BAND==KANE) kf = SMH[Material][IV]*sqrt(finalenergy*(1.+alphaK[Material][IV]*finalenergy));
       if(CONDUCTION_BAND==PARABOLIC) kf=SMH[Material][IV]*sqrt(finalenergy);
       cs = 1.-2.*rnd();
       sn = sqrt(1.-cs*cs);
       fai = 2.*PI*rnd();
       KX = kf*cs;
       KY = kf*sn*cos(fai);
       KZ = kf*sn*sin(fai);
       return;
      }
// Acoustic phonon
    if((r1<=SWK[Material][1][5][ie]) && j==0){
       finalenergy=superparticle_energy;
//       if(finalenergy<0.) return;
       finalk=sqrt(ksquared);
       j=1;
// determination of the final states
       if(CONDUCTION_BAND==KANE) kf = SMH[Material][IV]*sqrt(finalenergy*(1.+alphaK[Material][IV]*finalenergy));
       if(CONDUCTION_BAND==PARABOLIC) kf=SMH[Material][IV]*sqrt(finalenergy);
       cs = 1.-2.*rnd();
       sn = sqrt(1.-cs*cs);
       fai = 2.*PI*rnd();
       KX = kf*cs;
       KY = kf*sn*cos(fai);
       KZ = kf*sn*sin(fai);
       return;
    }
// Impurity scattering
    if((r1<=SWK[Material][1][6][ie]) && j==0){
     finalenergy=superparticle_energy;
//     if(finalenergy<0.) return;
     r2=rnd();
     cb=1.-r2/(0.5+(1.-r2)*ksquared/QD2);
     kf=ki;
// linea 30 -- determination of the final states
     sb=sqrt(1.-cb*cb);
     fai=2.*PI*rnd();
     cf=cos(fai);
     sf=sin(fai);
     skk=sqrt(KX*KX+KY*KY);
     a11=KY/skk;
     a12=KX*KZ/skk/ki;
     a13=KX/ki;
     a21=-KX/skk;
     a22=KY*KZ/skk/ki;
     a23=KY/ki;
     a32=-skk/ki;
     a33=KZ/ki;
     x1=kf*sb*cf;
     x2=kf*sb*sf;
     x3=kf*cb;
     KX=a11*x1+a12*x2+a13*x3;
     KY=a21*x1+a22*x2+a23*x3;
     KZ=a32*x2+a33*x3;
     return;
    }
    if((finalenergy<=0.) || j==0) return;
   }
// Selection of scattering process in the L-Valley
// ===================================================
  if(IV==2){
    r1 = rnd();
// Non-Polar optical phonons
// Emission of an optical phonon
      if(r1<=SWK[Material][2][1][ie] && j==0){
         finalenergy=superparticle_energy-HWO[Material][0];
         if(finalenergy<=0.) return;
         j=1;
// linea 20
       if(CONDUCTION_BAND==KANE) kf=SMH[Material][IV]*sqrt(finalenergy*(1.+alphaK[Material][IV]*finalenergy));
       if(CONDUCTION_BAND==PARABOLIC) kf=SMH[Material][IV]*sqrt(finalenergy);
       f=2.*ki*kf/(ki-kf)/(ki-kf);
       if(f<=0.) return;
       cb=(1.+f-pow((1.+2.*f),rnd()))/f;
// linea 30 -- determination of the final states
       sb=sqrt(1.-cb*cb);
       fai=2.*PI*rnd();
       cf=cos(fai);
       sf=sin(fai);
       skk=sqrt(KX*KX+KY*KY);
       a11=KY/skk;
       a12=KX*KZ/skk/ki;
       a13=KX/ki;
       a21=-KX/skk;
       a22=KY*KZ/skk/ki;
       a23=KY/ki;
       a32=-skk/ki;
       a33=KZ/ki;
       x1=kf*sb*cf;
       x2=kf*sb*sf;
       x3=kf*cb;
       KX=a11*x1+a12*x2+a13*x3;
       KY=a21*x1+a22*x2+a23*x3;
       KZ=a32*x2+a33*x3;
       return;
      }
// Absorbation of an optical phonon
      if((r1<=SWK[Material][2][2][ie]) && j==0){
         finalenergy=superparticle_energy+HWO[Material][0];
         j=1;
// linea 20
       if(CONDUCTION_BAND==KANE) kf=SMH[Material][IV]*sqrt(finalenergy*(1.+alphaK[Material][IV]*finalenergy));
       if(CONDUCTION_BAND==PARABOLIC) kf=SMH[Material][IV]*sqrt(finalenergy);
       f=2.*ki*kf/(ki-kf)/(ki-kf);
       if(f<=0.) return;
       cb=(1.+f-pow((1.+2.*f),rnd()))/f;
// linea 30 -- determination of the final states
       sb=sqrt(1.-cb*cb);
       fai=2.*PI*rnd();
       cf=cos(fai);
       sf=sin(fai);
       skk=sqrt(KX*KX+KY*KY);
       a11=KY/skk;
       a12=KX*KZ/skk/ki;
       a13=KX/ki;
       a21=-KX/skk;
       a22=KY*KZ/skk/ki;
       a23=KY/ki;
       a32=-skk/ki;
       a33=KZ/ki;
       x1=kf*sb*cf;
       x2=kf*sb*sf;
       x3=kf*cb;
       KX=a11*x1+a12*x2+a13*x3;
       KY=a21*x1+a22*x2+a23*x3;
       KZ=a32*x2+a33*x3;
       return;
      }
// Emission
    if((r1<=SWK[Material][2][3][ie]) && j==0){
       finalenergy=superparticle_energy-HWO[Material][0];
       if(finalenergy<=0.) return;
       j=1;
// determination of the final states
       if(CONDUCTION_BAND==KANE) kf = SMH[Material][IV]*sqrt(finalenergy*(1.+alphaK[Material][IV]*finalenergy));
       if(CONDUCTION_BAND==PARABOLIC) kf=SMH[Material][IV]*sqrt(finalenergy);
       cs = 1.-2.*rnd();
       sn = sqrt(1.-cs*cs);
       fai = 2.*PI*rnd();
       KX = kf*cs;
       KY = kf*sn*cos(fai);
       KZ = kf*sn*sin(fai);
       return;
    }
// Absorbation
    if((r1<=SWK[Material][2][4][ie]) && j==0){
       finalenergy=superparticle_energy+HWO[Material][0];
       if(finalenergy<=0.) return;
       j=1;
// determination of the final states
       if(CONDUCTION_BAND==KANE) kf = SMH[Material][IV]*sqrt(finalenergy*(1.+alphaK[Material][IV]*finalenergy));
       if(CONDUCTION_BAND==PARABOLIC) kf=SMH[Material][IV]*sqrt(finalenergy);
       cs = 1.-2.*rnd();
       sn = sqrt(1.-cs*cs);
       fai = 2.*PI*rnd();
       KX = kf*cs;
       KY = kf*sn*cos(fai);
       KZ = kf*sn*sin(fai);
       return;
    }
// Emission
    if((r1<=SWK[Material][2][5][ie]) && j==0){
       finalenergy=superparticle_energy-HWO[Material][0]+EMIN[Material][2]-EMIN[Material][1];
       if(finalenergy<=0.) return;
       IV=1;
       j=1;
// determination of the final states
       if(CONDUCTION_BAND==KANE) kf = SMH[Material][IV]*sqrt(finalenergy*(1.+alphaK[Material][IV]*finalenergy));
       if(CONDUCTION_BAND==PARABOLIC) kf=SMH[Material][IV]*sqrt(finalenergy);
       cs = 1.-2.*rnd();
       sn = sqrt(1.-cs*cs);
       fai = 2.*PI*rnd();
       KX = kf*cs;
       KY = kf*sn*cos(fai);
       KZ = kf*sn*sin(fai);
       return;
    }
// Absorbation
    if((r1<=SWK[Material][2][6][ie]) && j==0){
       finalenergy=superparticle_energy+HWO[Material][0]+EMIN[Material][2]-EMIN[Material][1];
       if(finalenergy<=0.) return;
       IV=1;
       j=1;
// determination of the final states
       if(CONDUCTION_BAND==KANE) kf = SMH[Material][IV]*sqrt(finalenergy*(1.+alphaK[Material][IV]*finalenergy));
       if(CONDUCTION_BAND==PARABOLIC) kf=SMH[Material][IV]*sqrt(finalenergy);
       cs = 1.-2.*rnd();
       sn = sqrt(1.-cs*cs);
       fai = 2.*PI*rnd();
       KX = kf*cs;
       KY = kf*sn*cos(fai);
       KZ = kf*sn*sin(fai);
       return;
    }
// Acoustic phonon emission
    if((r1<=SWK[Material][2][7][ie]) && j==0){
       finalenergy=superparticle_energy;
       finalk=sqrt(ksquared);
       j=1;
// determination of the final states
       if(CONDUCTION_BAND==KANE) kf = SMH[Material][IV]*sqrt(finalenergy*(1.+alphaK[Material][IV]*finalenergy));
       if(CONDUCTION_BAND==PARABOLIC) kf=SMH[Material][IV]*sqrt(finalenergy);
       cs = 1.-2.*rnd();
       sn = sqrt(1.-cs*cs);
       fai = 2.*PI*rnd();
       KX = kf*cs;
       KY = kf*sn*cos(fai);
       KZ = kf*sn*sin(fai);
       return;
    }
// Impurity scattering
    if((r1<=SWK[Material][2][8][ie]) && j==0){
     finalenergy=superparticle_energy;
     r2=rnd();
     cb=1.-r2/(0.5+(1.-r2)*ksquared/QD2);
     kf=ki;
// linea 30 -- determination of the final states
     sb=sqrt(1.-cb*cb);
     fai=2.*PI*rnd();
     cf=cos(fai);
     sf=sin(fai);
     skk=sqrt(KX*KX+KY*KY);
     a11=KY/skk;
     a12=KX*KZ/skk/ki;
     a13=KX/ki;
     a21=-KX/skk;
     a22=KY*KZ/skk/ki;
     a23=KY/ki;
     a32=-skk/ki;
     a33=KZ/ki;
     x1=kf*sb*cf;
     x2=kf*sb*sf;
     x3=kf*cb;
     KX=a11*x1+a12*x2+a13*x3;
     KY=a21*x1+a22*x2+a23*x3;
     KZ=a32*x2+a33*x3;
     return;     
    }
    if((finalenergy<=0.) || j==0) return;
   }
  }
}

// ===============================================
