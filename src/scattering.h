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

void scatter(particle_t *particle, int material)
{
    int has_scattered = 0,
        i  = 0,
        ie = 0;
    real ksquared = 0.,
         thesquareroot = 0.,
         superparticle_energy = 0.,
         r1 = 0.,
         finalenergy = 0.,
         finalk = 0.,
         cosinus = 0.,
         sinus = 0.,
         fai = 0.;
    real f   = 0.,
         cb  = 0.,
         cf  = 0.,
         sf  = 0.,
         skk = 0.,
         a11 = 0.,
         a12 = 0.,
         a13 = 0.,
         a21 = 0.,
         a22 = 0.,
         a23 = 0.,
         a32 = 0.,
         a33 = 0.,
         x1  = 0.,
         x2  = 0.,
         x3  = 0.,
         sb  = 0.,
         r2  = 0.,
         ki  = 0.,
         kf  = 0.,
         cs  = 0.,
         sn  = 0.;

    if(!mc_does_particle_exist(particle)) { return; }


    // ########################################
    // One-valley material Scattering Selection
    // ########################################
    if(NOVALLEY[material] == 1) {
        ksquared = mc_particle_ksquared(particle);

        if(g_config->conduction_band==FULL){
            real k, k2, k4;
            k = sqrt(ksquared) * 0.5 / PI * 1.e-12;
            k2 =  k * k;
            k4 = k2 * k2;
            k = sqrt(k2);
            // periodicity on reciprocal lattice
            superparticle_energy = CB_FULL[material][0] * k4 * k4 * k2
                                 + CB_FULL[material][1] * k4 * k4 * k
                                 + CB_FULL[material][2] * k4 * k4
                                 + CB_FULL[material][3] * k4 * k2 * k
                                 + CB_FULL[material][4] * k4 * k2
                                 + CB_FULL[material][5] * k4 * k
                                 + CB_FULL[material][6] * k4
                                 + CB_FULL[material][7] * k2 * k
                                 + CB_FULL[material][8] * k2
                                 + CB_FULL[material][9] * k
                                 + CB_FULL[material][10]; // in eV
        }
        if(g_config->conduction_band == KANE) {
            thesquareroot = sqrt(1. + 4. * alphaK[material][1] * HHM[material][0] * ksquared);
            superparticle_energy = (thesquareroot - 1.) / (2. * alphaK[material][1]);
        }
        if(g_config->conduction_band == PARABOLIC) {
            superparticle_energy = HHM[material][0] * ksquared; // in eV
        }

        if(superparticle_energy <= 0.) { return; }
        ie = ((int)(superparticle_energy / DE)) + 1;
        if(ie > DIME) { ie = DIME; }


        // ===============================
        // Selection of scattering process
        // ===============================
        r1 = rnd();

        // =========================
        // Non-Polar optical phonons
        for(i = 1; i <= 6; i++) {

            // Emission of an optical phonon
            if(r1 <= SWK[material][0][i*2-1][ie] && !has_scattered) {
                finalenergy = superparticle_energy - HWO[material][i-1];
                if(finalenergy <= 0.) { return; }
                has_scattered = 1;
            }

            // Absorbation of an optical phonon
            if((r1 <= SWK[material][0][i*2][ie]) && !has_scattered) {
                finalenergy = superparticle_energy + HWO[material][i-1];
                if(finalenergy <= 0.) { return; }
                has_scattered = 1;
            }
        }

        // =========================
        // Acoustic phonon
        if((r1 <= SWK[material][0][13][ie]) && !has_scattered) {
            finalenergy = superparticle_energy;
            if(finalenergy <= 0.) { return; }
            has_scattered = 1;
        }

        if((finalenergy <= 0.) || !has_scattered) { return; }


        // =================================
        // Determination of the final states
        // =================================
        if(g_config->conduction_band == FULL) {
            // look for a final k
            // bisection algorithm, probably not best but at least something to start from..
            real x = 0.0,
                 a = 0.0,
                 b = 1.0 / LATTCONST[material] * 1.e-12;
            real a2, a4, x4;
            real fx, fa;
            int k = 0;

            for(k = 0; k <= 512; k++) {
                x = 0.5 * (a + b);
                x2 = x * x;
                x4 = x2 * x2;
                a2 = a * a;
                a4 = a2 * a2;

                fa = (CB_FULL[material][0] * a4 * a4 * a2
                   +  CB_FULL[material][1] * a4 * a4 * a
                   +  CB_FULL[material][2] * a4 * a4
                   +  CB_FULL[material][3] * a4 * a2 * a
                   +  CB_FULL[material][4] * a4 * a2
                   +  CB_FULL[material][5] * a4 * a
                   +  CB_FULL[material][6] * a4
                   +  CB_FULL[material][7] * a2 * a
                   +  CB_FULL[material][8] * a2
                   +  CB_FULL[material][9] * a
                   +  CB_FULL[material][10])
                   -  finalenergy;
                fx = (CB_FULL[material][0] * x4 * x4 * x2
                   +  CB_FULL[material][1] * x4 * x4 * x
                   +  CB_FULL[material][2] * x4 * x4
                   +  CB_FULL[material][3] * x4 * x2 * x
                   +  CB_FULL[material][4] * x4 * x2
                   +  CB_FULL[material][5] * x4 * x
                   +  CB_FULL[material][6] * x4
                   +  CB_FULL[material][7] * x2 * x
                   +  CB_FULL[material][8] * x2
                   +  CB_FULL[material][9] * x
                   +  CB_FULL[material][10])
                   -  finalenergy;
                if((fa * fx) < 0.) { b = x; }
                else { a = x; }
            }
            finalk = x * 1.e12 * 2. * PI;
        }
        if(g_config->conduction_band == KANE) {
            finalk = SMH[material][0]
                   * sqrt(finalenergy * (1. + alphaK[material][1] * finalenergy));
        }
        if(g_config->conduction_band == PARABOLIC) {
            finalk = SMH[material][0] * sqrt(finalenergy);
        }

        cosinus = 1. - 2. * rnd();
        sinus = sqrt(1. - cosinus * cosinus);
        fai = 2. * PI * rnd();
        particle->kx = finalk * cosinus;
        particle->ky = finalk * sinus * cos(fai);
        particle->kz = finalk * sinus * sin(fai);
        return;
    }


    // ########################################
    // Two-valley material Scattering Selection
    // ########################################
    if(NOVALLEY[material] == 2) {
        ksquared = mc_particle_ksquared(particle);
        ki = sqrt(ksquared);

        if(g_config->conduction_band == KANE) {
            thesquareroot = sqrt(1. + 4. * alphaK[material][particle->valley]
                                         * HHM[material][particle->valley]
                                         * ksquared);
            superparticle_energy = (thesquareroot - 1.)
                                 / (2. * alphaK[material][particle->valley]);
        }
        if(g_config->conduction_band == PARABOLIC) {
            superparticle_energy = HHM[material][particle->valley] * ksquared;
        }

        if(superparticle_energy <= 0.) { return; }
        ie=((int)(superparticle_energy / DE)) + 1;
        if(ie > DIME) { ie=DIME; }


        // ===================================================
        // Selection of scattering process in the GAMMA-Valley
        // ===================================================
        if(particle->valley == 1) {
            r1 = rnd();

            // Non-Polar optical phonons
            // Emission of an optical phonon
            if(r1 <= SWK[material][1][1][ie] && !has_scattered) {
                finalenergy = superparticle_energy - HWO[material][0];
                if(finalenergy <= 0.) { return; }
                has_scattered = 1;

                if(g_config->conduction_band == KANE) {
                    kf = SMH[material][particle->valley]
                       * sqrt(finalenergy * (1. + alphaK[material][particle->valley] * finalenergy));
                }
                if(g_config->conduction_band == PARABOLIC) {
                    kf = SMH[material][particle->valley] * sqrt(finalenergy);
                }

                f = 2. * ki * kf / (ki - kf) / (ki - kf);
                if(f <= 0.) { return; }
                cb = (1. + f - pow(1. + 2. * f, rnd())) / f;

                // determination of the final states
                sb  = sqrt(1.-cb*cb);
                fai = 2.*PI*rnd();
                cf  = cos(fai);
                sf  = sin(fai);
                skk = sqrt(particle->kx*particle->kx
                         + particle->ky*particle->ky);
                a11 = particle->ky / skk;
                a12 = particle->kx * particle->kz / skk / ki;
                a13 = particle->kx / ki;
                a21 = -particle->kx / skk;
                a22 = particle->ky * particle->kz / skk / ki;
                a23 = particle->ky / ki;
                a32 = -skk / ki;
                a33 = particle->kz / ki;
                 x1 = kf * sb * cf;
                 x2 = kf * sb * sf;
                 x3 = kf * cb;
                particle->kx = a11 * x1 + a12 * x2 + a13 * x3;
                particle->ky = a21 * x1 + a22 * x2 + a23 * x3;
                particle->kz =            a32 * x2 + a33 * x3;
                return;
            }

            // Absorption of an optical phonon
            if((r1 <= SWK[material][1][2][ie]) && !has_scattered) {
                finalenergy = superparticle_energy + HWO[material][0];
                if(finalenergy <= 0.) { return; }
                has_scattered = 1;

                if(g_config->conduction_band == KANE) {
                    kf = SMH[material][particle->valley]
                       * sqrt(finalenergy * (1. + alphaK[material][particle->valley] * finalenergy));
                }
                if(g_config->conduction_band == PARABOLIC) {
                    kf=SMH[material][particle->valley] * sqrt(finalenergy);
                }

                f = 2. * ki * kf / (ki - kf) / (ki - kf);
                if(f <= 0.) { return; }
                cb = (1. + f - pow((1. + 2. * f), rnd())) / f;
                // determination of the final states
                sb=sqrt(1.-cb*cb);
                fai=2.*PI*rnd();
                cf=cos(fai);
                sf=sin(fai);
                skk=sqrt(particle->kx * particle->kx
                    +    particle->ky * particle->ky);
                a11 = particle->ky / skk;
                a12 = particle->kx * particle->kz / skk / ki;
                a13 = particle->kx / ki;
                a21 = -particle->kx / skk;
                a22 = particle->ky * particle->kz / skk / ki;
                a23 = particle->ky / ki;
                a32 = -skk / ki;
                a33 = particle->kz / ki;
                 x1 = kf * sb * cf;
                 x2 = kf * sb * sf;
                 x3 = kf * cb;
                particle->kx = a11 * x1 + a12 * x2 + a13 * x3;
                particle->ky = a21 * x1 + a22 * x2 + a23 * x3;
                particle->kz =            a32 * x2 + a33 * x3;
                return;
            }

            if((r1 <= SWK[material][1][3][ie]) && !has_scattered) {
                printf("Invalid Scatter!\n");
                return;
            }
            if((r1 <= SWK[material][1][4][ie]) && !has_scattered) {
                printf("Invalid Scatter!\n");
                return;
            }

            // Emission
            if((r1 <= SWK[material][1][5][ie]) && !has_scattered) {
                finalenergy = superparticle_energy - HWO[material][0]
                            + EMIN[material][1] - EMIN[material][2];
                if(finalenergy <= 0.) { return; }
                particle->valley = 2;
                has_scattered = 1;

                // determination of the final states
                if(g_config->conduction_band == KANE) {
                    kf = SMH[material][particle->valley]
                       * sqrt(finalenergy * (1. + alphaK[material][particle->valley] * finalenergy));
                }
                if(g_config->conduction_band == PARABOLIC) {
                    kf = SMH[material][particle->valley] * sqrt(finalenergy);
                }

                cs = 1. - 2. * rnd();
                sn = sqrt(1. - cs * cs);
                fai = 2. * PI * rnd();
                particle->kx = kf * cs;
                particle->ky = kf * sn * cos(fai);
                particle->kz = kf * sn * sin(fai);
                return;
            }

            // Absorption
            if((r1 <= SWK[material][1][6][ie]) && !has_scattered) {
                finalenergy = superparticle_energy + HWO[material][0]
                            + EMIN[material][1] - EMIN[material][2];
                if(finalenergy <= 0.) { return; }
                particle->valley = 2;
                has_scattered = 1;

                // determination of the final states
                if(g_config->conduction_band == KANE) {
                    kf = SMH[material][particle->valley]
                       * sqrt(finalenergy * (1. + alphaK[material][particle->valley] * finalenergy));
                }
                if(g_config->conduction_band == PARABOLIC) {
                    kf = SMH[material][particle->valley] * sqrt(finalenergy);
                }

                cs = 1. - 2. * rnd();
                sn = sqrt(1. - cs * cs);
                fai = 2. * PI * rnd();
                particle->kx = kf * cs;
                particle->ky = kf * sn * cos(fai);
                particle->kz = kf * sn * sin(fai);
                return;
            }

            // Acoustic phonon
            if((r1 <= SWK[material][1][7][ie]) && !has_scattered) {
                finalenergy = superparticle_energy;
                if(finalenergy <= 0.) { return; }
                finalk = sqrt(ksquared);
                has_scattered = 1;

                // determination of the final states
                if(g_config->conduction_band == KANE) {
                    kf = SMH[material][particle->valley]
                       * sqrt(finalenergy * (1. + alphaK[material][particle->valley] * finalenergy));
                }
                if(g_config->conduction_band == PARABOLIC) {
                    kf = SMH[material][particle->valley] * sqrt(finalenergy);
                }

                cs = 1. - 2. * rnd();
                sn = sqrt(1. - cs * cs);
                fai = 2. * PI * rnd();
                particle->kx = kf * cs;
                particle->ky = kf * sn * cos(fai);
                particle->kz = kf * sn * sin(fai);
                return;
            }

            // Impurity scattering
            if((r1 <= SWK[material][1][8][ie]) && !has_scattered) {
                finalenergy = superparticle_energy;
                if(finalenergy <= 0.) { return; }
                has_scattered = 1;

                r2 = rnd();
                cb = 1. - r2 / (0.5 + (1. - r2) * ksquared / QD2);
                kf = ki;

                // determination of the final states
                sb = sqrt(1. - cb * cb);
                fai = 2. * PI * rnd();
                cf = cos(fai);
                sf = sin(fai);
                skk = sqrt(particle->kx * particle->kx +
                           particle->ky * particle->ky);
                a11 = particle->ky / skk;
                a12 = particle->kx * particle->kz / skk / ki;
                a13 = particle->kx / ki;
                a21 = -particle->kx / skk;
                a22 = particle->ky * particle->kz / skk / ki;
                a23 = particle->ky / ki;
                a32 = -skk / ki;
                a33 = particle->kz / ki;
                x1 = kf * sb * cf;
                x2 = kf * sb * sf;
                x3 = kf * cb;
                particle->kx = a11 * x1 + a12 * x2 + a13 * x3;
                particle->ky = a21 * x1 + a22 * x2 + a23 * x3;
                particle->kz =            a32 * x2 + a33 * x3;
                return;
            }

            // Piezoelectric scattering
            if((r1 <= SWK[material][1][9][ie]) && !has_scattered) {
                finalenergy = superparticle_energy;
                if(finalenergy <= 0.) { return; }
                finalk = sqrt(ksquared);
                has_scattered = 1;

                // determination of the final states
                if(g_config->conduction_band == KANE) {
                    kf = SMH[material][particle->valley]
                       * sqrt(finalenergy * (1. + alphaK[material][particle->valley] * finalenergy));
                }
                if(g_config->conduction_band == PARABOLIC) {
                    kf = SMH[material][particle->valley] * sqrt(finalenergy);
                }

                cs = 1. - 2. * rnd();
                sn = sqrt(1. - cs * cs);
                fai = 2. * PI * rnd();
                particle->kx = kf * cs;
                particle->ky = kf * sn * cos(fai);
                particle->kz = kf * sn * sin(fai);
                return;
            }

            if((finalenergy <= 0.) || !has_scattered) { return; }
        }


        // ===============================================
        // Selection of scattering process in the L-Valley
        // ===============================================
        if(particle->valley == 2) {
            r1 = rnd();

            // Non-Polar optical phonons
            // Emission of an optical phonon
            if(r1 <= SWK[material][2][1][ie] && !has_scattered) {
                finalenergy = superparticle_energy - HWO[material][0];
                if(finalenergy <= 0.) { return; }
                has_scattered = 1;

                if(g_config->conduction_band == KANE) {
                    kf = SMH[material][particle->valley]
                       * sqrt(finalenergy * (1. + alphaK[material][particle->valley] * finalenergy));
                }
                if(g_config->conduction_band == PARABOLIC) {
                    kf = SMH[material][particle->valley] * sqrt(finalenergy);
                }

                f = 2. * ki * kf / (ki - kf) / (ki - kf);
                if(f <= 0.) { return; }
                cb = (1. + f - pow((1. + 2. * f), rnd())) / f;

                // determination of the final states
                sb = sqrt(1. - cb * cb);
                fai = 2. * PI * rnd();
                cf = cos(fai);
                sf = sin(fai);
                skk = sqrt(particle->kx * particle->kx +
                           particle->ky * particle->ky);
                a11 = particle->ky / skk;
                a12 = particle->kx * particle->kz / skk / ki;
                a13 = particle->kx / ki;
                a21 = -particle->kx / skk;
                a22 = particle->ky * particle->kz / skk / ki;
                a23 = particle->ky / ki;
                a32 = -skk / ki;
                a33 = particle->kz / ki;
                x1 = kf * sb * cf;
                x2 = kf * sb * sf;
                x3 = kf * cb;
                particle->kx = a11 * x1 + a12 * x2 + a13 * x3;
                particle->ky = a21 * x1 + a22 * x2 + a23 * x3;
                particle->kz =            a32 * x2 + a33 * x3;
                return;
            }

            // Absorbation of an optical phonon
            if((r1 <= SWK[material][2][2][ie]) && !has_scattered) {
                finalenergy = superparticle_energy + HWO[material][0];
                if(finalenergy <= 0.) { return; }
                has_scattered = 1;

                if(g_config->conduction_band == KANE) {
                    kf = SMH[material][particle->valley]
                       * sqrt(finalenergy * (1. + alphaK[material][particle->valley] * finalenergy));
                }
                if(g_config->conduction_band == PARABOLIC) {
                    kf = SMH[material][particle->valley] * sqrt(finalenergy);
                }

                f = 2. * ki * kf / (ki - kf) / (ki - kf);
                if(f <= 0.) { return; }
                cb = (1. + f - pow((1. + 2. * f), rnd())) / f;

                // determination of the final states
                sb = sqrt(1. - cb * cb);
                fai = 2. * PI * rnd();
                cf = cos(fai);
                sf = sin(fai);
                skk = sqrt(particle->kx * particle->kx +
                           particle->ky * particle->ky);
                a11 = particle->ky / skk;
                a12 = particle->kx * particle->kz / skk / ki;
                a13 = particle->kx / ki;
                a21 = -particle->kx / skk;
                a22 = particle->ky * particle->kz / skk / ki;
                a23 = particle->ky / ki;
                a32 = -skk / ki;
                a33 = particle->kz / ki;
                x1 = kf * sb * cf;
                x2 = kf * sb * sf;
                x3 = kf * cb;
                particle->kx = a11 * x1 + a12 * x2 + a13 * x3;
                particle->ky = a21 * x1 + a22 * x2 + a23 * x3;
                particle->kz =            a32 * x2 + a33 * x3;
                return;
            }


            // Emission
            if((r1 <= SWK[material][2][3][ie]) && !has_scattered) {
                finalenergy = superparticle_energy - HWO[material][0]
                            + EMIN[material][2] - EMIN[material][1];
                if(finalenergy <= 0.) { return; }
                particle->valley = 1;
                has_scattered = 1;

                // determination of the final states
                if(g_config->conduction_band == KANE) {
                    kf = SMH[material][particle->valley]
                       * sqrt(finalenergy * (1. + alphaK[material][particle->valley] * finalenergy));
                }
                if(g_config->conduction_band == PARABOLIC) {
                    kf = SMH[material][particle->valley] * sqrt(finalenergy);
                }

                cs = 1. - 2. * rnd();
                sn = sqrt(1. - cs * cs);
                fai = 2. * PI * rnd();
                particle->kx = kf * cs;
                particle->ky = kf * sn * cos(fai);
                particle->kz = kf * sn * sin(fai);
                return;
            }

            // Absorbation
            if((r1 <= SWK[material][2][4][ie]) && !has_scattered) {
                finalenergy = superparticle_energy + HWO[material][0]
                            + EMIN[material][2] - EMIN[material][1];
                if(finalenergy <= 0.) { return; }
                particle->valley = 1;
                has_scattered = 1;

                // determination of the final states
                if(g_config->conduction_band == KANE) {
                    kf = SMH[material][particle->valley]
                       * sqrt(finalenergy * (1. + alphaK[material][particle->valley] * finalenergy));
                }
                if(g_config->conduction_band == PARABOLIC) {
                    kf = SMH[material][particle->valley] * sqrt(finalenergy);
                }

                cs = 1. - 2. * rnd();
                sn = sqrt(1. - cs * cs);
                fai = 2. * PI * rnd();
                particle->kx = kf * cs;
                particle->ky = kf * sn * cos(fai);
                particle->kz = kf * sn * sin(fai);
                return;
            }

            // Emission
            if((r1 <= SWK[material][2][5][ie]) && !has_scattered) {
                finalenergy = superparticle_energy - HWO[material][0];
                if(finalenergy <= 0.) { return; }
                has_scattered = 1;

                if(g_config->conduction_band == KANE) {
                    kf = SMH[material][particle->valley]
                       * sqrt(finalenergy * (1. + alphaK[material][particle->valley] * finalenergy));
                }
                if(g_config->conduction_band == PARABOLIC) {
                    kf = SMH[material][particle->valley] * sqrt(finalenergy);
                }

                cs = 1. - 2. * rnd();
                sn = sqrt(1. - cs * cs);
                fai = 2. * PI * rnd();
                particle->kx = kf * cs;
                particle->ky = kf * sn * cos(fai);
                particle->kz = kf * sn * sin(fai);
                return;
            }

            // Absorbation
            if((r1 <= SWK[material][2][6][ie]) && !has_scattered) {
                finalenergy = superparticle_energy + HWO[material][0];
                if(finalenergy <= 0.) { return; }
                has_scattered = 1;

                // determination of the final states
                if(g_config->conduction_band==KANE) {
                    kf = SMH[material][particle->valley]
                       * sqrt(finalenergy * (1. + alphaK[material][particle->valley] * finalenergy));
                }
                if(g_config->conduction_band==PARABOLIC) {
                    kf = SMH[material][particle->valley] * sqrt(finalenergy);
                }

                cs = 1. - 2. * rnd();
                sn = sqrt(1. - cs * cs);
                fai = 2. * PI * rnd();
                particle->kx = kf * cs;
                particle->ky = kf * sn * cos(fai);
                particle->kz = kf * sn * sin(fai);
                return;
            }

            // Acoustic phonon emission
            if((r1 <= SWK[material][2][7][ie]) && !has_scattered) {
                finalenergy = superparticle_energy;
                if(finalenergy <= 0.) { return; }
                has_scattered = 1;

                finalk = sqrt(ksquared);
                // determination of the final states
                if(g_config->conduction_band == KANE) {
                    kf = SMH[material][particle->valley]
                       * sqrt(finalenergy * (1. + alphaK[material][particle->valley] * finalenergy));
                }
                if(g_config->conduction_band == PARABOLIC) {
                    kf = SMH[material][particle->valley] * sqrt(finalenergy);
                }
                cs = 1. - 2. * rnd();
                sn = sqrt(1. - cs * cs);
                fai = 2. * PI * rnd();
                particle->kx = kf * cs;
                particle->ky = kf * sn * cos(fai);
                particle->kz = kf * sn * sin(fai);
                return;
            }

            // Impurity scattering
            if((r1 <= SWK[material][2][8][ie]) && !has_scattered) {
                finalenergy = superparticle_energy;
                if(finalenergy <= 0.) { return; }
                has_scattered = 1;

                r2 = rnd();
                cb = 1. - r2 / (0.5 + (1. - r2) * ksquared / QD2);
                kf = ki;

                // determination of the final states
                sb = sqrt(1. - cb * cb);
                fai=2. * PI * rnd();
                cf = cos(fai);
                sf = sin(fai);
                skk = sqrt(particle->kx * particle->kx +
                           particle->ky * particle->ky);
                a11 = particle->ky / skk;
                a12 = particle->kx * particle->kz / skk / ki;
                a13 = particle->kx / ki;
                a21 = -particle->kx / skk;
                a22 = particle->ky * particle->kz / skk / ki;
                a23 = particle->ky / ki;
                a32 = -skk / ki;
                a33 = particle->kz / ki;
                x1 = kf * sb * cf;
                x2 = kf * sb * sf;
                x3 = kf * cb;
                particle->kx = a11 * x1 + a12 * x2 + a13 * x3;
                particle->ky = a21 * x1 + a22 * x2 + a23 * x3;
                particle->kz =            a32 * x2 + a33 * x3;
                return;
            }

            // Piezoelectric scattering
            if((r1 <= SWK[material][2][9][ie]) && !has_scattered) {
                finalenergy = superparticle_energy;
                if(finalenergy <= 0.) { return; }
                finalk = sqrt(ksquared);
                has_scattered = 1;

                // determination of the final states
                if(g_config->conduction_band == KANE) {
                    kf = SMH[material][particle->valley]
                       * sqrt(finalenergy * (1. + alphaK[material][particle->valley] * finalenergy));
                }
                if(g_config->conduction_band == PARABOLIC) {
                    kf = SMH[material][particle->valley] * sqrt(finalenergy);
                }

                cs = 1. - 2. * rnd();
                sn = sqrt(1. - cs * cs);
                fai = 2. * PI * rnd();
                particle->kx = kf * cs;
                particle->ky = kf * sn * cos(fai);
                particle->kz = kf * sn * sin(fai);
                return;
            }

            if((finalenergy <= 0.) || !has_scattered) { return; }
        }
    }
}
