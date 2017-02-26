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

int scatter(Particle *particle, Material *material) {
    int has_scattered = 0;
    double ksquared = 0.,
           ki = 0.,
           kf = 0.,
           superparticle_energy = 0.,
           finalenergy = 0.;

    if(!mc_does_particle_exist(particle)) { return has_scattered; }


    // ########################################
    // One-valley material Scattering Selection
    // ########################################
    if(material->cb.num_valleys == 1) {
        ksquared = mc_particle_ksquared(particle);
        ki = sqrt(ksquared);

        superparticle_energy = mc_particle_energy(particle);

        if(superparticle_energy <= 0.) { return has_scattered; }
        int ie = ((int)(superparticle_energy / DE)) + 1;
        if(ie > DIME) { ie = DIME; }


        // ===============================
        // Selection of scattering process
        // ===============================
        double r1 = rnd();

        // =========================
        // Non-Polar optical phonons
        for(int i = 1; i <= 6; i++) {

            // Emission of an optical phonon
            if(r1 <= SWK[material->id][0][i*2-1][ie] && !has_scattered) {
                finalenergy = superparticle_energy - material->hwo[i-1];
                if(finalenergy <= 0.) { return has_scattered; }
                has_scattered = 1;
            }

            // Absorption of an optical phonon
            if((r1 <= SWK[material->id][0][i*2][ie]) && !has_scattered) {
                finalenergy = superparticle_energy + material->hwo[i-1];
                if(finalenergy <= 0.) { return has_scattered; }
                has_scattered = 1;
            }
        }

        // =========================
        // Acoustic phonon
        if((r1 <= SWK[material->id][0][13][ie]) && !has_scattered) {
            finalenergy = superparticle_energy;
            if(finalenergy <= 0.) { return has_scattered; }
            has_scattered = 1;
        }

        if((finalenergy <= 0.) || !has_scattered) { return has_scattered; }


        // =================================
        // Determination of the final states
        // =================================
        mc_calculate_isotropic_k(particle, finalenergy);
        return has_scattered;
    }


    // ########################################
    // Two-valley material Scattering Selection
    // ########################################
    if(material->cb.num_valleys == 2) {
        ksquared = mc_particle_ksquared(particle);
        ki = sqrt(ksquared);

        superparticle_energy = mc_particle_energy(particle);

        if(superparticle_energy <= 0.) { return has_scattered; }
        int ie = ((int)(superparticle_energy / DE)) + 1;
        if(ie > DIME) { ie = DIME; }


        // ===================================================
        // Selection of scattering process in the GAMMA-Valley
        // ===================================================
        if(particle->valley == 1) {
            double r1 = rnd();

            // NPOP Emission
            if(r1 <= SWK[material->id][1][1][ie] && !has_scattered) {
                finalenergy = superparticle_energy - material->hwo[0];
                if(finalenergy <= 0.) { return has_scattered; }
                has_scattered = 1;

                if(g_config->conduction_band == KANE) {
                    kf = material->cb.smh[particle->valley]
                       * sqrt(finalenergy * (1. + material->cb.alpha[particle->valley] * finalenergy));
                }
                if(g_config->conduction_band == PARABOLIC) {
                    kf = material->cb.smh[particle->valley] * sqrt(finalenergy);
                }

                double f = 2. * ki * kf / (ki - kf) / (ki - kf);
                if(f <= 0.) { return has_scattered; }
                double cb = (1. + f - pow(1. + 2. * f, rnd())) / f;

                mc_calculate_anisotropic_k(particle, ki, kf, cb);
                return has_scattered;
            }

            // NPOP Absorption
            if((r1 <= SWK[material->id][1][2][ie]) && !has_scattered) {
                finalenergy = superparticle_energy + material->hwo[0];
                if(finalenergy <= 0.) { return has_scattered; }
                has_scattered = 1;

                if(g_config->conduction_band == KANE) {
                    kf = material->cb.smh[particle->valley]
                       * sqrt(finalenergy * (1. + material->cb.alpha[particle->valley] * finalenergy));
                }
                if(g_config->conduction_band == PARABOLIC) {
                    kf=material->cb.smh[particle->valley] * sqrt(finalenergy);
                }

                double f = 2. * ki * kf / (ki - kf) / (ki - kf);
                if(f <= 0.) { return has_scattered; }
                double cb = (1. + f - pow((1. + 2. * f), rnd())) / f;

                mc_calculate_anisotropic_k(particle, ki, kf, cb);
                return has_scattered;
            }

            // POP Emission 1->1 (not allowed)
            if((r1 <= SWK[material->id][1][3][ie]) && !has_scattered) {
                printf("Invalid Scatter!\n");
                return has_scattered;
            }
            // POP Absorption 1->1 (not allowed)
            if((r1 <= SWK[material->id][1][4][ie]) && !has_scattered) {
                printf("Invalid Scatter!\n");
                return has_scattered;
            }

            // POP Emission 1->2
            if((r1 <= SWK[material->id][1][5][ie]) && !has_scattered) {
                finalenergy = superparticle_energy - material->hwo[0]
                            + (material->cb.emin[1] - material->cb.emin[2]);
                if(finalenergy <= 0.) { return has_scattered; }
                particle->valley = 2;
                has_scattered = 1;

                mc_calculate_isotropic_k(particle, finalenergy);
                return has_scattered;
            }

            // POP Absorption 1->2
            if((r1 <= SWK[material->id][1][6][ie]) && !has_scattered) {
                finalenergy = superparticle_energy + material->hwo[0]
                            + (material->cb.emin[1] - material->cb.emin[2]);
                if(finalenergy <= 0.) { return has_scattered; }
                particle->valley = 2;
                has_scattered = 1;

                mc_calculate_isotropic_k(particle, finalenergy);
                return has_scattered;
            }

            // Acoustic phonon
            if((r1 <= SWK[material->id][1][7][ie]) && !has_scattered) {
                finalenergy = superparticle_energy;
                if(finalenergy <= 0.) { return has_scattered; }
                kf = sqrt(ksquared);
                has_scattered = 1;

                mc_calculate_isotropic_k(particle, finalenergy);
                return has_scattered;
            }

            // Impurity scattering
            if((r1 <= SWK[material->id][1][8][ie]) && !has_scattered) {
                finalenergy = superparticle_energy;
                if(finalenergy <= 0.) { return has_scattered; }
                has_scattered = 1;

                double r2 = rnd();
                double cb = 1. - r2 / (0.5 + (1. - r2) * ksquared / QD2);
                kf = ki;

                mc_calculate_anisotropic_k(particle, ki, kf, cb);
                return has_scattered;
            }

            // Piezoelectric scattering
            if((r1 <= SWK[material->id][1][9][ie]) && !has_scattered) {
                finalenergy = superparticle_energy;
                if(finalenergy <= 0.) { return has_scattered; }
                kf = sqrt(ksquared);
                has_scattered = 1;

                mc_calculate_isotropic_k(particle, finalenergy);
                return has_scattered;
            }

            if((finalenergy <= 0.) || !has_scattered) { return has_scattered; }
        }


        // ===============================================
        // Selection of scattering process in the L-Valley
        // ===============================================
        if(particle->valley == 2) {
            double r1 = rnd();

            // NPOP Emission
            if(r1 <= SWK[material->id][2][1][ie] && !has_scattered) {
                finalenergy = superparticle_energy - material->hwo[0];
                if(finalenergy <= 0.) { return has_scattered; }
                has_scattered = 1;

                if(g_config->conduction_band == KANE) {
                    kf = material->cb.smh[particle->valley]
                       * sqrt(finalenergy * (1. + material->cb.alpha[particle->valley] * finalenergy));
                }
                if(g_config->conduction_band == PARABOLIC) {
                    kf = material->cb.smh[particle->valley] * sqrt(finalenergy);
                }

                double f = 2. * ki * kf / (ki - kf) / (ki - kf);
                if(f <= 0.) { return has_scattered; }
                double cb = (1. + f - pow((1. + 2. * f), rnd())) / f;

                mc_calculate_anisotropic_k(particle, ki, kf, cb);
                return has_scattered;
            }

            // NPOP Absorption
            if((r1 <= SWK[material->id][2][2][ie]) && !has_scattered) {
                finalenergy = superparticle_energy + material->hwo[0];
                if(finalenergy <= 0.) { return has_scattered; }
                has_scattered = 1;

                if(g_config->conduction_band == KANE) {
                    kf = material->cb.smh[particle->valley]
                       * sqrt(finalenergy * (1. + material->cb.alpha[particle->valley] * finalenergy));
                }
                if(g_config->conduction_band == PARABOLIC) {
                    kf = material->cb.smh[particle->valley] * sqrt(finalenergy);
                }

                double f = 2. * ki * kf / (ki - kf) / (ki - kf);
                if(f <= 0.) { return has_scattered; }
                double cb = (1. + f - pow((1. + 2. * f), rnd())) / f;

                mc_calculate_anisotropic_k(particle, ki, kf, cb);
                return has_scattered;
            }


            // POP Emission 2->1
            if((r1 <= SWK[material->id][2][3][ie]) && !has_scattered) {
                finalenergy = superparticle_energy - material->hwo[0]
                            + (material->cb.emin[2] - material->cb.emin[1]);
                if(finalenergy <= 0.) { return has_scattered; }
                particle->valley = 1;
                has_scattered = 1;

                mc_calculate_isotropic_k(particle, finalenergy);
                return has_scattered;
            }

            // POP Absorption 2->1
            if((r1 <= SWK[material->id][2][4][ie]) && !has_scattered) {
                finalenergy = superparticle_energy + material->hwo[0]
                            + (material->cb.emin[2] - material->cb.emin[1]);
                if(finalenergy <= 0.) { return has_scattered; }
                particle->valley = 1;
                has_scattered = 1;

                mc_calculate_isotropic_k(particle, finalenergy);
                return has_scattered;
            }

            // POP Emission 2->2
            if((r1 <= SWK[material->id][2][5][ie]) && !has_scattered) {
                finalenergy = superparticle_energy - material->hwo[0];
                if(finalenergy <= 0.) { return has_scattered; }
                has_scattered = 1;

                mc_calculate_isotropic_k(particle, finalenergy);
                return has_scattered;
            }

            // POP Absorption 2->2
            if((r1 <= SWK[material->id][2][6][ie]) && !has_scattered) {
                finalenergy = superparticle_energy + material->hwo[0];
                if(finalenergy <= 0.) { return has_scattered; }
                has_scattered = 1;

                mc_calculate_isotropic_k(particle, finalenergy);
                return has_scattered;
            }

            // Acoustic phonon
            if((r1 <= SWK[material->id][2][7][ie]) && !has_scattered) {
                finalenergy = superparticle_energy;
                if(finalenergy <= 0.) { return has_scattered; }
                has_scattered = 1;

                mc_calculate_isotropic_k(particle, finalenergy);
                return has_scattered;
            }

            // Impurity scattering
            if((r1 <= SWK[material->id][2][8][ie]) && !has_scattered) {
                finalenergy = superparticle_energy;
                if(finalenergy <= 0.) { return has_scattered; }
                has_scattered = 1;

                double r2 = rnd();
                double cb = 1. - r2 / (0.5 + (1. - r2) * ksquared / QD2);
                kf = ki;

                mc_calculate_anisotropic_k(particle, ki, kf, cb);
                return has_scattered;
            }

            // Piezoelectric scattering
            if((r1 <= SWK[material->id][2][9][ie]) && !has_scattered) {
                finalenergy = superparticle_energy;
                if(finalenergy <= 0.) { return has_scattered; }
                kf = sqrt(ksquared);
                has_scattered = 1;

                mc_calculate_isotropic_k(particle, finalenergy);
                return has_scattered;
            }

            if((finalenergy <= 0.) || !has_scattered) { return has_scattered; }
        }
    }

    return has_scattered;
}
