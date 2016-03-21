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

void MCparameters(int material)
{
    real wo,no,dos,aco,oge[7],oga[7];
    real dos1,dos2,am1,am2,cl,deq,dij;
    real hwe,hwij,wij,we,ne,nij;
    real poe,poa,ope,opa,eqe,eqa,qmin,qmax;
    real initialenergy,sei,finalenergy,sef;
    real eps,epf,ep,bimp,cimp,qd;
    real ak,qq,wk;
    int ie,i;
    real z2=4., zf;
    real gamma1_initial, gamma1_final,
         gamma2_initial, gamma2_final,
         gamma_initial,  gamma_final,
         sqgamma_initial, sqgamma_final;
    real overlapA, overlapB, overlapC, overlap, rate;

    ISEED = 38467.;  //  initial value for random number generator

    // These definitions are valid for every material
    BKTQ=KB*TL/Q; // in eV
    QH=Q/HBAR;

    // Material with 2 valleys
    // #######################
    if(NOVALLEY[material] == 2) {
        // Effective mass for the GAMMA and L-Valley, respectively
        am1 = MSTAR[material][1] * M;
        am2 = MSTAR[material][2] * M;
        // Dielectric constant - static, hi-freq and combination
        eps = EPSR[material] * EPS0;
        epf = EPF[material]  * EPS0;
        ep  = 1. / ((1. / epf) - (1. / eps));


        // =======================
        // == Phonon Scattering ==
        // =======================
        cl = RHO[material] * pow(UL[material], 2.);  // rho * vs^2
        deq = dij = DTK[material][0] * Q; // deformation potential
        hwe = hwij = HWO[material][0]; // phonon energy

        SMH[material][1] = sqrt(2. * am1 * Q) / HBAR;
        SMH[material][2] = sqrt(2. * am2 * Q) / HBAR;
        HHM[material][1] = HBAR * HBAR / (2. * am1 * Q);
        HHM[material][2] = HBAR * HBAR / (2. * am2 * Q);
        HM[material][1]  = HBAR / am1;
        HM[material][2]  = HBAR / am2;

        // Phonon frequency
        wo  = HWO[material][0] * Q / HBAR;
        wij = hwij * Q / HBAR;
        we  = hwe  * Q / HBAR;

        // Phonon number
        no  = 1. / (exp(HWO[material][0] / BKTQ) - 1.);
        nij = 1. / (exp(hwij / BKTQ) - 1.);
        ne  = 1. / (exp(hwe  / BKTQ) - 1.);

        // Effective mass and density of states prefactor for each valley
        real mstar[4];
        real dos[4];
        real alpha[4];
        // band structure parameters
        for(int v = 1; v <= NOVALLEY[material]; v++) {
            mstar[v] = MSTAR[material][v] * M;
            dos[v] = pow(sqrt(2. * mstar[v]) / HBAR, 3.) / (2. * PI * PI);
            if(CONDUCTION_BAND == PARABOLIC) { alpha[v] = 0.; }
            else if(CONDUCTION_BAND == KANE) { alpha[v] = alphaK[material][v]; }
        }


        poe = Q * Q * wo / (4. * PI * ep) * (no + 1.);
        poa = poe * no / (no + 1.);

        // Polar Optical Prefactor
        //   Emission
        poe= Q / 8. / PI / ep * Q * wo * ( no + 1.);
        //   Absorption
        poa= poe * no / (1. + no);

        // Acoustic Prefactor (Elastic)
        aco = 2. * PI * DA[material] / Q * DA[material] * BKTQ / HBAR * Q / cl;

        // Non-polar Optical Prefactor
        //   Emission
        ope = PI * dij / wij * dij / RHO[material] / Q * (nij + 1.);
        //   Absorption
        opa = ope * nij / (1. + nij);
        //   Emission
        eqe = PI * deq / we * deq / RHO[material] / Q * (ne + 1.);
        //   Absorption
        eqa= eqe * ne / (1. + ne);

        // =========================
        // == Impurity Scattering ==
        // =========================
        cimp = CIMP; // impurity concentration
        QD2 = Q * cimp / BKTQ / eps;
        qd = sqrt(QD2);
        bimp = 2. * PI * cimp * Q * Q / HBAR * Q / eps / eps;


        // =====================================
        // == Calculation of scattering rates ==
        // =====================================
        char filename[150];
        FILE *pop_emission,  *pop_absorption,
             *npop_emission, *npop_absorption,
             *ac_elastic,    *imp_elastic;
        char f_band_model[20];
        char f_material[20];

        if(CONDUCTION_BAND == PARABOLIC) { sprintf(f_band_model, "para"); }
        else if(CONDUCTION_BAND == KANE) { sprintf(f_band_model, "kane"); }
        if(material == GAAS) { sprintf(f_material, "GaAs"); }
        else if(material == SILICON) { sprintf(f_material, "Si"); }
        else if(material == GAN) { sprintf(f_material, "GaN"); }

        sprintf(filename, "pop_emission_%s-%s.csv", f_material, f_band_model);
        pop_emission = fopen(filename, "w");
        fprintf(pop_emission, "energy,rate\n");

        sprintf(filename, "pop_absorption_%s-%s.csv", f_material, f_band_model);
        pop_absorption = fopen(filename, "w");
        fprintf(pop_absorption, "energy,rate\n");

        sprintf(filename, "npop_emission_%s-%s.csv", f_material, f_band_model);
        npop_emission = fopen(filename, "w");
        fprintf(npop_emission, "energy,rate\n");

        sprintf(filename, "npop_absorption_%s-%s.csv", f_material, f_band_model);
        npop_absorption = fopen(filename, "w");
        fprintf(npop_absorption, "energy,rate\n");

        sprintf(filename, "ac_elastic_%s-%s.csv", f_material, f_band_model);
        ac_elastic = fopen(filename, "w");
        fprintf(ac_elastic, "energy,rate\n");

        sprintf(filename, "imp_elastic_%s-%s.csv", f_material, f_band_model);
        imp_elastic = fopen(filename, "w");
        fprintf(imp_elastic, "energy,rate\n");

        for(ie = 1; ie <= DIME; ie++) {
            initialenergy = DE * (real)(ie);
            sei = sqrt(initialenergy);

            for(int v = 1; v < NOVALLEY[material]; ++v) {

                // ===============================
                // == Optical Phonon Scattering ==
                // ===============================
                if(OPTICALPHONONS == ON) {

                    SWK[material][v][0][ie] = 0.;

                    // Polar Optical Phonons
                    // ---------------------
                    for(int i = 1; i <= 2; ++i) {
                        SWK[material][v][i][ie] = SWK[material][v][i-1][ie];

                        real prefactor = 0.;
                        FILE *f;
                        // Emission
                        if(i == 1) {
                            finalenergy = initialenergy - HWO[material][0];
                            prefactor = poe;
                            f = pop_emission;
                        }
                        // Absorption
                        if(i == 2) {
                            finalenergy = initialenergy + HWO[material][0];
                            prefactor = poa;
                            f = pop_absorption;
                        }

                        if(finalenergy <= 0.) { // Negative energy, ignore
                            fprintf(f, "%g,%g\n", initialenergy, 0.);
                            SWK[material][1][i][ie] += 0.;
                        }
                        else {
                            // non-parabolicity precalculations
                            //   γ  = E(1 + αE) = E*γ1
                            //   γ1 = 1 + αE
                            //   γ2 = 1 + 2αE
                            gamma1_initial = 1. +      alpha[v] * initialenergy;
                            gamma1_final   = 1. +      alpha[v] * finalenergy;
                            gamma2_initial = 1. + 2. * alpha[v] * initialenergy;
                            gamma2_final   = 1. + 2. * alpha[v] * finalenergy;
                            gamma_initial = initialenergy * gamma1_initial;
                            gamma_final   = finalenergy   * gamma1_final;
                            sqgamma_initial = sqrt(gamma_initial);
                            sqgamma_final   = sqrt(gamma_final);
                            sef = sqrt(finalenergy);
                            qmax = sqgamma_initial + sqgamma_final;
                            qmin = fabs(sqgamma_initial - sqgamma_final);
                            overlapA = pow(2. * gamma1_initial * gamma1_final
                                           + alpha[v] * (gamma_initial + gamma_final), 2.);
                            overlapB = 2. * alpha[v] * sqgamma_initial * sqgamma_final
                                       * (4. * gamma1_initial * gamma1_final
                                          + alpha[v] * (gamma_initial + gamma_final));
                            overlapC = 4. * gamma1_initial * gamma1_final
                                          * gamma2_initial * gamma2_final;
                            overlap = (overlapA * log(qmax / qmin) - overlapB) / overlapC;

                            rate = prefactor * SMH[material][v] * gamma2_final / sqgamma_initial * overlap / Q;
                            fprintf(f, "%g,%g\n", initialenergy, rate);
                            SWK[material][v][i][ie] += rate;
                        }
                    } // POP scattering

                    // Non-polar Optical Phonons - Intervalley
                    // ---------------------------------------
                    for(int v2 = 1; v2 <= NOVALLEY[material]; ++v2) {
                        // scatter from valley v -> v2
                        // when v == v2, scattering to equivalent valley
                        for(int i = 3; i <= 4; ++i) {
                            SWK[material][v][i][ie] = SWK[material][v][i-1][ie];

                            real prefactor = 0.;
                            FILE *f;

                            // Emission
                            if(i == 3) {
                                finalenergy = initialenergy - hwij - fabs(EMIN[material][v] - EMIN[material][v2]);
                                prefactor = ope * Q;
                                f = npop_emission;
                            }
                            // Absorption
                            if(i == 4) {
                                finalenergy = initialenergy + hwij - fabs(EMIN[material][v] - EMIN[material][v2]);
                                prefactor = opa * Q;
                                f = npop_absorption;
                            }

                            if(finalenergy <= 0.) { // Negative energy, ignore
                                if(v != v2) {
                                    fprintf(f, "%g,%g\n", initialenergy, 0.);
                                }
                                SWK[material][v][i][ie] += 0.;
                            }
                            else {
                                gamma1_initial = 1. +      alpha[v]  * initialenergy;
                                gamma1_final   = 1. +      alpha[v2] * finalenergy;
                                gamma2_initial = 1. + 2. * alpha[v]  * initialenergy;
                                gamma2_final   = 1. + 2. * alpha[v2] * finalenergy;
                                gamma_initial = Q * initialenergy * gamma1_initial;
                                gamma_final   = Q * finalenergy   * gamma1_final;
                                sqgamma_initial = sqrt(gamma_initial);
                                sqgamma_final   = sqrt(gamma_final);
                                overlap = (gamma1_initial * gamma1_final) / (gamma2_initial * gamma2_final);
                                zf = ZSCATTER[material][v][v2];

                                rate = prefactor * zf * dos[v2] * sqgamma_final * gamma2_final * overlap;
                                if(v != v2) {
                                    fprintf(f, "%g,%g\n", initialenergy, rate);
                                }
                                SWK[material][v][i][ie] += rate;
                            }
                        } // NPOP scattering
                    } // secondary loop over valleys
                }
                else { // no optical phonon scattering
                    for(int i = 1; i <= 4; ++i) {
                        SWK[material][v][i][ie] = 0.;
                    }
                } // optical phonon scattering

                // ================================
                // == Acoustic Phonon Scattering ==
                // ================================
                if(ACOUSTICPHONONS == ON) {
                    finalenergy = initialenergy; // elastic scattering assumption

                    real gamma1 = 1. +      alpha[v] * finalenergy;
                    real gamma2 = 1. + 2. * alpha[v] * finalenergy;
                    real gamma = finalenergy * Q * gamma1;
                    real sqgamma = sqrt(gamma);

                    overlapA = gamma1 * gamma1;
                    overlapB = (alpha[v] * alpha[v] * finalenergy * finalenergy) / 3;
                    overlapC = gamma2 * gamma2;
                    overlap = (overlapA + overlapB) / overlapC;

                    rate = aco * Q * dos[v] * sqgamma * gamma2 * overlap;

                    fprintf(ac_elastic, "%g,%g\n", initialenergy, rate);
                    SWK[material][v][5][ie] = SWK[material][v][4][ie] + rate;
                }
                else { // no acoustic phonon scattering
                    SWK[material][v][5][ie] = 0.;
                    fprintf(ac_elastic, "%g,%g\n", initialenergy, 0.);
                } // acoustic phonon scattering

                // =========================
                // == Impurity Scattering ==
                // =========================
                if(IMPURITYPHONONS == ON) {
                    finalenergy = initialenergy; // elastic scattering

                    real gamma1 = 1. +      alpha[v] * finalenergy;
                    real gamma2 = 1. + 2. * alpha[v] * finalenergy;
                    real gamma = finalenergy * Q * gamma1;
                    real sqgamma = sqrt(gamma);

                    real prefactor = 2 * PI * CIMP * Q * Q * Q * Q / ( HBAR * eps * eps);
                    real qd2 = Q * Q * CIMP / (BKTQ * Q * eps);
                    real k2 = 2 * mstar[v] * gamma / (HBAR * HBAR);
                    real screening = qd2 * (4. * k2 + qd2);

                    rate = prefactor * sqgamma * gamma2 * dos[v] / screening;
                    fprintf(imp_elastic, "%g,%g\n", initialenergy, rate);
                    SWK[material][v][6][ie] = SWK[material][v][5][ie] + rate;
                }
                else { // no impurity scattering
                    SWK[material][v][6][ie] = 0.;
                    fprintf(imp_elastic, "%g,%g\n", initialenergy, 0.);
                } // impurity scattering

            } // loop over valleys


// Impurity scattering
   if(IMPURITYPHONONS==ON){
    finalenergy=initialenergy;
    sef=sqrt(finalenergy*(1.+alphaK[material][1]*finalenergy));
    ak=SMH[material][1]*sef;
    qq=QD2*(4.*ak*ak+QD2);
    wk=bimp/qq*sef*dos1*(1.+2.*alphaK[material][1]*finalenergy);
    if(wk>1.e14) wk=1.e14;
    SWK[material][1][6][ie]=SWK[material][1][5][ie]+wk;
   }
   else {
    // NO IMPURITY SCATTERING
    SWK[material][1][6][ie]=SWK[material][1][5][ie]+0.0;
   }
// L-valley
// ========
   if(OPTICALPHONONS==ON){
// Polar optical phonon
// Emission
    finalenergy=initialenergy-HWO[material][0];
    if(finalenergy>0.){
      sef=sqrt(finalenergy);
      qmax=sef+sei;
      qmin=sei-sef;
      SWK[material][2][1][ie]=poe*SMH[material][2]*sei/initialenergy/Q*log(qmax/qmin);
    }
    else SWK[material][2][1][ie]=0.;
// Absorption
    finalenergy=initialenergy+HWO[material][0];
    sef=sqrt(finalenergy);
    qmax=sef+sei;
    qmin=sef-sei;
    SWK[material][2][2][ie]=SWK[material][2][1][ie]
                 +poa*SMH[material][2]*sei/initialenergy/Q*log(qmax/qmin);
// Non-polar optical phonon
// Emission
    finalenergy=initialenergy-hwe;
    if(finalenergy>0.){
      sef=sqrt(finalenergy*(1.+alphaK[material][2]*finalenergy));
      SWK[material][2][3][ie]=SWK[material][2][2][ie]
                   +(z2-1.)*eqe*sef*dos2*(1.+2.*alphaK[material][2]*finalenergy);
    }
    else SWK[material][2][3][ie]=SWK[material][2][2][ie];
// Absorption
    finalenergy=initialenergy+hwe;
    sef=sqrt(finalenergy*(1.+alphaK[material][2]*finalenergy));
    SWK[material][2][4][ie]=SWK[material][2][3][ie]
                 +(z2-1.)*eqa*sef*dos2*(1.+2.*alphaK[material][2]*finalenergy);

// Emission
    finalenergy=initialenergy-hwij+EMIN[material][2]-EMIN[material][1];
    if(finalenergy>0.){
     sef=sqrt(finalenergy*(1.+alphaK[material][2]*finalenergy));
     SWK[material][2][5][ie]=SWK[material][2][4][ie]
                            +ope*sef*dos1*(1.+2.*alphaK[material][2]*finalenergy);
    }
    else SWK[material][2][5][ie]=SWK[material][2][4][ie];
// Absorption
    finalenergy=initialenergy+hwij+EMIN[material][2]-EMIN[material][1];
    if(finalenergy>0.){
     sef=sqrt(finalenergy*(1.+alphaK[material][2]*finalenergy));
     SWK[material][2][6][ie]=SWK[material][2][5][ie]
                            +opa*sef*dos1*(1.+2.*alphaK[material][2]*finalenergy);
    }
    else SWK[material][2][6][ie]=SWK[material][2][5][ie];
   }
   else {
      // NO OPTICAL SCATTERING
      SWK[material][2][1][ie]=0.0;
      SWK[material][2][2][ie]=0.0;
      SWK[material][2][3][ie]=0.0;
      SWK[material][2][4][ie]=0.0;
      SWK[material][2][5][ie]=0.0;
      SWK[material][2][6][ie]=0.0;
   }
   if(ACOUSTICPHONONS==ON){
// Acoustic phonon
    finalenergy=initialenergy;
    sef=sqrt(finalenergy*(1.+alphaK[material][2]*finalenergy));
    SWK[material][2][7][ie]=SWK[material][2][6][ie]
                           +aco*sef*dos2*(1.+2.*alphaK[material][2]*finalenergy);
   }
   else {
      // NO ACOUSTIC SCATTERING
      SWK[material][2][7][ie]=SWK[material][2][6][ie]+0.0;
   }
// ==================
// --->  } <--- This is a bug! Removed by J.M.Sellier on 24 dec.2006 (SR)
// (Thanks to Kun-Yuan Xu who helped me to find it!)
// ==================
   if(IMPURITYPHONONS==ON){
// Impurity scattering
    finalenergy=initialenergy;
    sef=sqrt(finalenergy*(1.+alphaK[material][2]*finalenergy));
    ak=SMH[material][2]*sef;
    qq=QD2*(4.*ak*ak+QD2);
    wk=bimp/qq*sef*dos2*(1.+2.*alphaK[material][2]*finalenergy);
    if(wk>1.e14) wk=1.e14;
    SWK[material][2][8][ie]=SWK[material][2][7][ie]+wk;
   }
   else {
    // NO IMPURITY SCATTERING
    SWK[material][2][8][ie]=SWK[material][2][7][ie]+0.0;
   }
  } // <--- this is the correct place for the symbol } (J.M.Sellier,24/12/2006)
// Evalutation of gamma
  GM[material]=SWK[material][1][6][1];
  for(ie=1;ie<=DIME;ie++){
    if(SWK[material][1][6][ie]>GM[material]) GM[material]=SWK[material][1][6][ie];
    if(SWK[material][2][8][ie]>GM[material]) GM[material]=SWK[material][2][8][ie];
  }
  printf("GAMMA[%d] = %g \n",material, GM[material]);
  for(i=1;i<=6;i++)
    for(ie=1;ie<=DIME;ie++)
      SWK[material][1][i][ie]/=GM[material];
  for(i=1;i<=8;i++)
    for(ie=1;ie<=DIME;ie++)
      SWK[material][2][i][ie]/=GM[material];
 }
// End of two-valleys materials
// ############################

// Material = one valley
// #####################
 if(NOVALLEY[material]==1){
  SMH[material][0]=sqrt(2.*MSTAR[material][1]*M*Q)/HBAR;
  HHM[material][0]=HBAR*HBAR/(2.*MSTAR[material][1]*M*Q);
  HM[material][0]=HBAR/(MSTAR[material][1]*M);
// Density of states
  dos=pow((sqrt(2.*MSTAR[material][1]*M)*sqrt(Q)/HBAR),3.)/(4.*PI*PI);
// constant for the acoustic phonon
  aco=2.*PI*(DA[material]/Q)*DA[material]*(BKTQ/HBAR)
     *(Q/(RHO[material]*UL[material]*UL[material]));
// Constants for the 6 Silicon-like non-polar optical phonons
   for(i=1;i<=6;i++){
// i-th Optical Phonon
    oge[i]=0.;
    oga[i]=0.;
    if(ZF[material][i-1]!=0.){
      wo=HWO[material][i-1]*Q/HBAR; // frequency of phonon
      no=1./(exp(HWO[material][i-1]/BKTQ) - 1.); // population of phonons
      oge[i]=ZF[material][i-1]*PI*(DTK[material][i-1]*Q/wo)
            *((DTK[material][i-1]*Q/RHO[material])/Q)*(no+1.);
      oga[i]=oge[i]*no/(1.+no);
    }
   }
// Calculation of scattering rates
   for(ie=1; ie<=DIME; ++ie) SWK[material][0][0][ie]=0.;
   for(ie=1; ie<=DIME; ++ie){
    initialenergy=DE*((real) ie);
    sei=sqrt(initialenergy);
    if(OPTICALPHONONS==ON){
// non polar optical phonons
     for(i=1;i<=6;i++){
      finalenergy=initialenergy-HWO[material][i-1];
      SWK[material][0][i*2-1][ie]=SWK[material][0][i*2-2][ie];
      if(finalenergy>0.){
       sef=sqrt(finalenergy*(1.+alphaK[material][1]*finalenergy));
       SWK[material][0][i*2-1][ie]=SWK[material][0][i*2-2][ie]
                     +oge[i]*sef*dos*(1.+2.*alphaK[material][1]*finalenergy);
      }
      finalenergy=initialenergy+HWO[material][i-1];
      SWK[material][0][i*2][ie]=SWK[material][0][i*2-1][ie];
      if(finalenergy>0.){
       sef=sqrt(finalenergy*(1.+alphaK[material][1]*finalenergy));
       SWK[material][0][i*2][ie]=SWK[material][0][i*2-1][ie]
                   +oga[i]*sef*dos*(1.+2.*alphaK[material][1]*finalenergy);
      }
     }
    }
    else {
     // NO OPTICAL SCATTERING
     SWK[material][0][1][ie]=0.0;
     SWK[material][0][2][ie]=0.0;
     SWK[material][0][3][ie]=0.0;
     SWK[material][0][4][ie]=0.0;
     SWK[material][0][5][ie]=0.0;
     SWK[material][0][6][ie]=0.0;
     SWK[material][0][7][ie]=0.0;
     SWK[material][0][8][ie]=0.0;
     SWK[material][0][9][ie]=0.0;
     SWK[material][0][10][ie]=0.0;
     SWK[material][0][11][ie]=0.0;
     SWK[material][0][12][ie]=0.0;
    }
    if(ACOUSTICPHONONS==ON){
// Acoustic phonon
    finalenergy=initialenergy;
    sef=sqrt(finalenergy*(1.+alphaK[material][1]*finalenergy));
    SWK[material][0][13][ie]=SWK[material][0][12][ie]
                  +aco*sef*dos*(1.+2.*alphaK[material][1]*finalenergy);
    }
    else {
    // NO ACOUSTIC PHONONS
    SWK[material][0][13][ie]=SWK[material][0][12][ie]+0.0;
    }
   }
// Evaluation of gamma
  GM[material]=SWK[material][0][13][1];
  for(ie=1;ie<=DIME;++ie)
    if(SWK[material][0][13][ie]>GM[material]) GM[material]=SWK[material][0][13][ie];
  printf("GAMMA[%d] = %g \n",material, GM[material]);
  for(ie=1;ie<=DIME;ie++)
    for(i=1;i<=13;i++)
      SWK[material][0][i][ie]/=GM[material];
 }
// End of one-valley material
}

// =======================================================
