/* scattering_rates.h -- This file is part of Archimedes release 1.1.0.
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
// Last modif. : 04 may.2016, J. Marini
// ######################################################

// Definition of the variables needed for taking into
// account scatterings from acoustic and optical,non-polar
// phonons (which are the most relevant scatterings in common semiconductors like Silicon)

int calculate_scattering_rates(Material *material) {
    real wo,no,aco,oge[7],oga[7];
    real cl,dij;
    real hwij,wij,nij;
    real poe,poa,ope,opa,qmin,qmax;
    real initialenergy,finalenergy,sef;
    real eps,epf,ep,cimp;
    int ie,i;
    real zf = 0.;
    real gamma1_initial, gamma1_final,
         gamma2_initial, gamma2_final,
         gamma_initial,  gamma_final,
         sqgamma_initial, sqgamma_final;
    real overlapA, overlapB, overlapC, overlap, rate;

    // These definitions are valid for every material
    BKTQ=KB*g_config->lattice_temp/Q; // in eV
    int num_valleys = material->cb.num_valleys;

    // Material with 2 valleys
    // #######################
    if(num_valleys >= 2) {
        // Dielectric constant - static, hi-freq and combination
        eps = material->eps_static * EPS0;
        epf = material->eps_hf  * EPS0;
        ep  = 1. / ((1. / epf) - (1. / eps));


        // =======================
        // == Phonon Scattering ==
        // =======================
        cl = material->rho * pow(material->ul, 2.);  // rho * vs^2
        dij = material->dtk[0] * Q; // deformation potential
        hwij = material->hwo[0]; // phonon energy

        // Phonon frequency
        wo  = material->hwo[0] * Q / HBAR;
        wij = hwij * Q / HBAR;

        // Phonon number
        no  = 1. / (exp(material->hwo[0] / BKTQ) - 1.);
        nij = 1. / (exp(hwij / BKTQ) - 1.);

        // Effective mass and density of states prefactor for each valley
        real dos[MAX_VALLEYS];
        // band structure parameters
        for(int v = 1; v <= num_valleys; v++) {
            dos[v] = pow(sqrt(2. * material->cb.mstar[v] * M) / HBAR, 3.) / (4. * PI * PI);
        }


        poe = Q * Q * wo / (4. * PI * ep) * (no + 1.);
        poa = poe * no / (no + 1.);

        // Polar Optical Prefactor
        //   Emission
        poe= Q / 8. / PI / ep * Q * wo * ( no + 1.);
        //   Absorption
        poa= poe * no / (1. + no);

        // Acoustic Prefactor (Elastic)
        aco = 2. * PI * material->da / Q * material->da * BKTQ / HBAR * Q / cl;

        // Non-polar Optical Prefactor
        //   Emission
        ope = PI * dij / wij * dij / material->rho / Q * (nij + 1.);
        //   Absorption
        opa = ope * nij / (1. + nij);

        // =========================
        // == Impurity Scattering ==
        // =========================
        cimp = g_config->impurity_conc; // impurity concentration
        QD2 = Q * cimp / BKTQ / eps;


        // =====================================
        // == Calculation of scattering rates ==
        // =====================================
        // scattering rate files, indexed by scattering type and valley number
        //   1: Polar Optical Emission
        //   2: Polar Optical Absorption
        //   3: Non-polar Optical Emission    V->1
        //   4: Non-polar Optical Absorption  V->1
        //   5: Non-polar Optical Emission    V->2
        //   6: Non-polar Optical Absorption  V->2
        //   7: Acoustic Phonon
        //   8: Impurity Scattering
        //   9: Piezoelectric Scattering

        char *f_band_model = mc_band_model_name(g_config->conduction_band);
        char *f_material = mc_material_name(material);

        int imax = 5 + 2 * num_valleys;

        FILE * scattering_rates[12][MAX_VALLEYS];
        if(g_config->scattering_output) {
            char *f_scattering[] = {
                "neutral_imp",
                "charged_imp",
                "ac_elastic",
                "piezo_elastic",
                "pop_emission",  "pop_absorption",
                "npop_emission", "npop_absorption",
                "npop_emission", "npop_absorption",
                "npop_emission", "npop_absorption"
            };
            char filename[150];

            for(int s = 0; s <= imax; ++s ) {
                for(int v = 1; v <= num_valleys; ++v) {
                    if(s >= 6) {
                        int v2 = (int)((s - 4) / 2);
                        sprintf(filename, "%s_%d-%d-%s-%s.csv", f_scattering[s], v, v2, f_material, f_band_model);
                        scattering_rates[s][v] = fopen(filename, "w");
                    }
                    else {
                        sprintf(filename, "%s_%d-%s-%s.csv", f_scattering[s], v, f_material, f_band_model);
                        scattering_rates[s][v] = fopen(filename, "w");
                    }
                }
            }
        }

        for(ie = 1; ie <= DIME; ie++) {
            initialenergy = DE * (real)(ie);

            for(int v = 1; v <= num_valleys; ++v) {
                SWK[material->id][v][0][ie] = 0.;

                // =========================
                // == Impurity Scattering ==
                // =========================
                if(g_config->neutral_impurity_scattering == ON) {
                    int i = 0;
                    finalenergy = initialenergy; // elastic scattering

                    double bohr = 4. * PI * eps * HBAR * HBAR / (material->vb.mstar[0] * M * Q * Q);
                    double a = 12.5 + HBAR*HBAR / (2. * material->cb.mstar[v] * M * KB * g_config->lattice_temp * bohr * bohr);
                    double c = 3.4 * g_config->neutral_impurity_conc * HBAR / (material->cb.mstar[v] * M) * pow(HBAR*HBAR / (2 * material->cb.mstar[v] * M * KB), 3. / 2.);

                    rate = c * pow(g_config->lattice_temp * a, -3. / 2.) / (bohr * bohr);
                    SWK[material->id][v][i][ie] = 0. + rate;
                    if(g_config->scattering_output) {
                        fprintf(scattering_rates[i][v], "%g,%g\n", initialenergy, rate);
                    }
                }
                else {
                    int i = 0;
                    SWK[material->id][v][i][ie] = 0.;
                    if(g_config->scattering_output) {
                        fprintf(scattering_rates[i][v], "%g,%g\n", initialenergy, 0.);
                    }
                }

                if(g_config->impurity_scattering == ON) {
                    int i = 1;
                    finalenergy = initialenergy; // elastic scattering

                    real gamma1 = 1. +      material->cb.alpha[v] * finalenergy;
                    real gamma2 = 1. + 2. * material->cb.alpha[v] * finalenergy;
                    real gamma = finalenergy * Q * gamma1;
                    real sqgamma = sqrt(gamma);

                    real prefactor = 2 * PI * g_config->impurity_conc * Q * Q * Q * Q / ( HBAR * eps * eps);
                    real qd2 = QD2;
                    if(g_config->thomas_fermi_screening == ON) {
                        real Nlh = g_config->impurity_conc / (1. + pow(material->vb.mstar[0] / material->vb.mstar[1], 1.5));
                        qd2 = Q * Q * material->vb.mstar[1] * M * pow(3.0 * Nlh, 0.33) / (pow(PI, 1.33) * eps * HBAR * HBAR);
                    }
                    real k2 = 2 * material->cb.mstar[v] * M * gamma / (HBAR * HBAR);
                    real screening = qd2 * (4. * k2 + qd2);

                    rate = prefactor * sqgamma * gamma2 * dos[v] / screening;
                    SWK[material->id][v][i][ie] = SWK[material->id][v][i-1][ie] + rate;
                    if(g_config->scattering_output) {
                        fprintf(scattering_rates[i][v], "%g,%g\n", initialenergy, rate);
                    }
                }
                else { // no impurity scattering
                    int i = 1;
                    SWK[material->id][v][i][ie] = SWK[material->id][v][i-1][ie];
                    if(g_config->scattering_output) {
                        fprintf(scattering_rates[i][v], "%g,%g\n", initialenergy, 0.);
                    }
                }

                // ================================
                // == Acoustic Phonon Scattering ==
                // ================================
                if(g_config->acoustic_phonon_scattering == ON) {
                    int i = 2;
                    finalenergy = initialenergy; // elastic scattering assumption

                    real gamma1 = 1. +      material->cb.alpha[v] * finalenergy;
                    real gamma2 = 1. + 2. * material->cb.alpha[v] * finalenergy;
                    real gamma = finalenergy * Q * gamma1;
                    real sqgamma = sqrt(gamma);

                    overlapA = gamma1 * gamma1;
                    overlapB = (material->cb.alpha[v] * material->cb.alpha[v] * finalenergy * finalenergy) / 3;
                    overlapC = gamma2 * gamma2;
                    overlap = (overlapA + overlapB) / overlapC;

                    rate = aco * Q * dos[v] * sqgamma * gamma2 * overlap;

                    SWK[material->id][v][i][ie] = SWK[material->id][v][i-1][ie] + rate;
                    if(g_config->scattering_output) {
                        fprintf(scattering_rates[i][v], "%g,%g\n", initialenergy, rate);
                    }
                }
                else { // no acoustic phonon scattering
                    int i = 2;
                    SWK[material->id][v][i][ie] = SWK[material->id][v][i-1][ie];
                    if(g_config->scattering_output) {
                        fprintf(scattering_rates[i][v], "%g,%g\n", initialenergy, 0.);
                    }
                }


                if(g_config->piezoelectric_scattering == ON) {
                    int i = 3;
                    finalenergy = initialenergy; // elastic scattering

                    real gamma1 = 1. +      material->cb.alpha[v] * finalenergy;
                    real gamma2 = 1. + 2. * material->cb.alpha[v] * finalenergy;
                    real gamma = finalenergy * Q * gamma1;
                    real sqgamma = sqrt(gamma);

                    real prefactor = Q*Q * material->kav * BKTQ * Q / (HBAR * HBAR * eps * 4 * PI);
                    real qd2 = QD2;
                    if(g_config->thomas_fermi_screening == ON) {
                        real Nlh = g_config->impurity_conc / (1. + pow(material->vb.mstar[0] / material->vb.mstar[1], 1.5));
                        qd2 = Q * Q * material->vb.mstar[1] * M * pow(3.0 * Nlh, 0.33) / (pow(PI, 1.33) * eps * HBAR * HBAR);
                    }
                    real a = material->cb.mstar[v] * M * finalenergy * Q / (HBAR * HBAR * qd2);
                    real screening = log(1. + a) - a / (1. + a);
                    real rate = prefactor * sqrt(material->cb.mstar[v] * M / 2.) * gamma2 * screening / sqgamma;
                    SWK[material->id][v][i][ie] = SWK[material->id][v][i-1][ie] + rate;
                    if(g_config->scattering_output) {
                        fprintf(scattering_rates[i][v], "%g,%g\n", initialenergy, rate);
                    }
                }
                else {
                    int i = 3;
                    SWK[material->id][v][i][ie] = SWK[material->id][v][i-1][ie];
                    if(g_config->scattering_output) {
                        fprintf(scattering_rates[i][v], "%g,%g\n", initialenergy, 0.);
                    }
                }

                // ===============================
                // == Optical Phonon Scattering ==
                // ===============================
                if(g_config->optical_phonon_scattering == ON) {
                    SWK[material->id][v][0][ie] = 0.;

                    // Polar Optical Phonons
                    // ---------------------
                    for(int i = 4; i <= 5; ++i) {
                        SWK[material->id][v][i][ie] = SWK[material->id][v][i-1][ie];

                        real prefactor = 0.;
                        // Emission
                        if(i == 4) {
                            finalenergy = initialenergy - material->hwo[0];
                            prefactor = poe;
                        }
                        // Absorption
                        if(i == 5) {
                            finalenergy = initialenergy + material->hwo[0];
                            prefactor = poa;
                        }

                        if(finalenergy <= 0.) { // Negative energy, ignore
                            SWK[material->id][v][i][ie] += 0.;
                            if(g_config->scattering_output) {
                                fprintf(scattering_rates[i][v], "%g,%g\n", initialenergy, 0.);
                            }
                        }
                        else {
                            // non-parabolicity precalculations
                            //   γ  = E(1 + αE) = E*γ1
                            //   γ1 = 1 + αE
                            //   γ2 = 1 + 2αE
                            gamma1_initial = 1. +      material->cb.alpha[v] * initialenergy;
                            gamma1_final   = 1. +      material->cb.alpha[v] * finalenergy;
                            gamma2_initial = 1. + 2. * material->cb.alpha[v] * initialenergy;
                            gamma2_final   = 1. + 2. * material->cb.alpha[v] * finalenergy;
                            gamma_initial = initialenergy * gamma1_initial;
                            gamma_final   = finalenergy   * gamma1_final;
                            sqgamma_initial = sqrt(gamma_initial);
                            sqgamma_final   = sqrt(gamma_final);
                            sef = sqrt(finalenergy);
                            qmax = sqgamma_initial + sqgamma_final;
                            qmin = fabs(sqgamma_initial - sqgamma_final);
                            overlapA = pow(2. * gamma1_initial * gamma1_final
                                           + material->cb.alpha[v] * (gamma_initial + gamma_final), 2.);
                            overlapB = 2. * material->cb.alpha[v] * sqgamma_initial * sqgamma_final
                                       * (4. * gamma1_initial * gamma1_final
                                          + material->cb.alpha[v] * (gamma_initial + gamma_final));
                            overlapC = 4. * gamma1_initial * gamma1_final
                                          * gamma2_initial * gamma2_final;
                            overlap = (overlapA * log(qmax / qmin) - overlapB) / overlapC;

                            rate = prefactor * material->cb.smh[v] * gamma2_final / sqgamma_initial * overlap / Q;
                            SWK[material->id][v][i][ie] += rate;
                            if(g_config->scattering_output) {
                                fprintf(scattering_rates[i][v], "%g,%g\n", initialenergy, rate);
                            }
                        }
                    } // POP scattering

                    // Non-polar Optical Phonons - Intervalley
                    // ---------------------------------------
                    for(int v2 = 1; v2 <= num_valleys; ++v2) {
                        // scatter from valley v -> v2
                        // when v == v2, scattering to equivalent valley
                        for(int n = 0; n <= 1; ++n) {
                            int i = 2 * v2 + n + 4;
                            SWK[material->id][v][i][ie] = SWK[material->id][v][i-1][ie];

                            real prefactor = 0.;

                            // Emission
                            if(n == 0) {
                                finalenergy = initialenergy - hwij + (material->cb.emin[v] - material->cb.emin[v2]);
                                prefactor = ope * Q;
                            }
                            // Absorption
                            if(n == 1) {
                                finalenergy = initialenergy + hwij + (material->cb.emin[v] - material->cb.emin[v2]);
                                prefactor = opa * Q;
                            }

                            if(finalenergy <= 0.) { // Negative energy, ignore
                                SWK[material->id][v][i][ie] += 0.;
                                if(g_config->scattering_output) {
                                    fprintf(scattering_rates[i][v], "%g,%g\n", initialenergy, 0.);
                                }
                            }
                            else {
                                gamma1_initial = 1. +      material->cb.alpha[v]  * initialenergy;
                                gamma1_final   = 1. +      material->cb.alpha[v2] * finalenergy;
                                gamma2_initial = 1. + 2. * material->cb.alpha[v]  * initialenergy;
                                gamma2_final   = 1. + 2. * material->cb.alpha[v2] * finalenergy;
                                gamma_initial = Q * initialenergy * gamma1_initial;
                                gamma_final   = Q * finalenergy   * gamma1_final;
                                sqgamma_initial = sqrt(gamma_initial);
                                sqgamma_final   = sqrt(gamma_final);
                                overlap = (gamma1_initial * gamma1_final) / (gamma2_initial * gamma2_final);
                                zf = material->zscatter[v][v2];

                                rate = prefactor * zf * dos[v2] * sqgamma_final * gamma2_final * overlap;
                                SWK[material->id][v][i][ie] += rate;
                                if(g_config->scattering_output) {
                                    fprintf(scattering_rates[i][v], "%g,%g\n", initialenergy, rate);
                                }
                            }
                        } // NPOP scattering
                    } // secondary loop over valleys
                }
                else { // no optical phonon scattering
                    for(int i = 6; i < 6 + 2 * num_valleys; ++i) {
                        SWK[material->id][v][i][ie] = SWK[material->id][v][i-1][ie];
                    }
                } // optical phonon scattering


            } // loop over valleys
        } // loop over energy


        // Evalutation of gamma
        GM[material->id] = 0.0;
        for(ie = 1; ie <= DIME; ie++) {
            for(int v = 1; v <= num_valleys; ++v) {
                if(SWK[material->id][v][imax][ie] > GM[material->id]) {
                    GM[material->id] = SWK[material->id][v][imax][ie];
                }
            }
        }
        printf("GAMMA[%s] = %g\n", f_material, GM[material->id]);
        for(ie = 1; ie <= DIME; ie++) {
            for(int v = 1; v <= num_valleys; ++v) {
                for(i = 0; i <= imax; i++) {
                    SWK[material->id][v][i][ie] /= GM[material->id];
                }
            }
        }

        if(g_config->scattering_output) {
            for(int i = 0; i <= imax; ++i) {
                for(int v = 1; v <= num_valleys; ++v) {
                    fclose(scattering_rates[i][v]);
                }
            }
        }
    }

    // End of two-valleys materials
    // ############################

// Material = one valley
// #####################
 if(num_valleys==1){
        char *f_material = mc_material_name(material);
// Density of states
  real dos=pow((sqrt(2.*material->cb.mstar[1]*M)*sqrt(Q)/HBAR),3.)/(4.*PI*PI);
// constant for the acoustic phonon
  aco=2.*PI*(material->da/Q)*material->da*(BKTQ/HBAR)
     *(Q/(material->rho*material->ul*material->ul));
// Constants for the 6 Silicon-like non-polar optical phonons
   for(i=1;i<=6;i++){
// i-th Optical Phonon
    oge[i]=0.;
    oga[i]=0.;
    if(material->zf[i-1]!=0.){
      wo=material->hwo[i-1]*Q/HBAR; // frequency of phonon
      no=1./(exp(material->hwo[i-1]/BKTQ) - 1.); // population of phonons
      oge[i]=material->zf[i-1]*PI*(material->dtk[i-1]*Q/wo)
            *((material->dtk[i-1]*Q/material->rho)/Q)*(no+1.);
      oga[i]=oge[i]*no/(1.+no);
    }
   }
// Calculation of scattering rates
   for(ie=1; ie<=DIME; ++ie) SWK[material->id][0][0][ie]=0.;
   for(ie=1; ie<=DIME; ++ie){
    initialenergy=DE*((real) ie);
    if(g_config->optical_phonon_scattering==ON){
// non polar optical phonons
     for(i=1;i<=6;i++){
      finalenergy=initialenergy-material->hwo[i-1];
      SWK[material->id][0][i*2-1][ie]=SWK[material->id][0][i*2-2][ie];
      if(finalenergy>0.){
       sef=sqrt(finalenergy*(1.+material->cb.alpha[1]*finalenergy));
       SWK[material->id][0][i*2-1][ie]=SWK[material->id][0][i*2-2][ie]
                     +oge[i]*sef*dos*(1.+2.*material->cb.alpha[1]*finalenergy);
      }
      finalenergy=initialenergy+material->hwo[i-1];
      SWK[material->id][0][i*2][ie]=SWK[material->id][0][i*2-1][ie];
      if(finalenergy>0.){
       sef=sqrt(finalenergy*(1.+material->cb.alpha[1]*finalenergy));
       SWK[material->id][0][i*2][ie]=SWK[material->id][0][i*2-1][ie]
                   +oga[i]*sef*dos*(1.+2.*material->cb.alpha[1]*finalenergy);
      }
     }
    }
    else {
     // NO OPTICAL SCATTERING
     SWK[material->id][0][1][ie]=0.0;
     SWK[material->id][0][2][ie]=0.0;
     SWK[material->id][0][3][ie]=0.0;
     SWK[material->id][0][4][ie]=0.0;
     SWK[material->id][0][5][ie]=0.0;
     SWK[material->id][0][6][ie]=0.0;
     SWK[material->id][0][7][ie]=0.0;
     SWK[material->id][0][8][ie]=0.0;
     SWK[material->id][0][9][ie]=0.0;
     SWK[material->id][0][10][ie]=0.0;
     SWK[material->id][0][11][ie]=0.0;
     SWK[material->id][0][12][ie]=0.0;
    }
    if(g_config->acoustic_phonon_scattering==ON){
// Acoustic phonon
    finalenergy=initialenergy;
    sef=sqrt(finalenergy*(1.+material->cb.alpha[1]*finalenergy));
    SWK[material->id][0][13][ie]=SWK[material->id][0][12][ie]
                  +aco*sef*dos*(1.+2.*material->cb.alpha[1]*finalenergy);
    }
    else {
    // NO ACOUSTIC PHONONS
    SWK[material->id][0][13][ie]=SWK[material->id][0][12][ie]+0.0;
    }
   }
// Evaluation of gamma
  GM[material->id]=SWK[material->id][0][13][1];
  for(ie=1;ie<=DIME;++ie)
    if(SWK[material->id][0][13][ie]>GM[material->id]) GM[material->id]=SWK[material->id][0][13][ie];
  printf("GAMMA[%s] = %g\n", f_material, GM[material->id]);
  for(ie=1;ie<=DIME;ie++)
    for(i=1;i<=13;i++)
      SWK[material->id][0][i][ie]/=GM[material->id];
 }
// End of one-valley material

    return 0;
}
