/* optical_absorption.h -- This file is part of Archimedes release 1.2.0.
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


// Returns the optical joint density of states for the given conduction and valence bands
// evaluated at the given energy
double optical_joint_DOS(int material, int conduction_band, int valence_band, double energy) {
    double mc = MSTAR[material][conduction_band],
           mv = MSTAR_VB[material][valence_band];
    double mr = mc * mv / (mc + mv);
    double alpha = alphaK[material][conduction_band];

    double deltaE = DELTAE_VB[material][valence_band] + EG[material];
    double e2 = (mr / mc) * (energy - deltaE);

    double gamma = Q * e2 * (1. + alpha * e2);
    double gamma2 = (1. + 2. * alpha * e2);

    double dos = mr * M * sqrt(2. * mr * M) * sqrt(gamma) * gamma2 / (PI * PI * HBAR * HBAR * HBAR);
    return dos;
}


// Returns the transition rate (1/s) between the given conduction and valence bands at the given
// energy
double optical_transition_rate(int material, int conduction_band,
                             int valence_band, double photon_energy) {
    if(photon_energy < g_materials[material].Eg + g_materials[material].cb[conduction_band].emin + g_materials[material].vb[valence_band].emin) {
        return 0.; // photon energy too small
    }

    double prefactor = PI * HBAR * Q * Q * EG[material]
                     / (3. * EPSR[material] * EPS0 * M * photon_energy);
    double rate = prefactor * optical_joint_DOS(material, conduction_band, valence_band, photon_energy);
    return rate;
}


// Returns the absorption coefficient at the given energy. Takes into account direct interband
// transitions only from HH, LH, SO valence bands to G-1 conduction band
double absorption_coefficient(int material, double photon_energy) {
    double rate = 0.;
    for(int v = 0; v < 3; ++v) {
        rate += optical_transition_rate(material, 1, v, photon_energy) * sqrt(EPSR[material]) / VLIGHT;
    }

    return rate * g_materials[material].abs_correction;
}


// Calculates the number of photoexcited electrons in the given node at the specified energy.
// Photons are assumed to enter device at x=0
// TODO: number of carriers needs to depend on energy - relative W @ E vs max W
int electrons_in_cell(Node *node, double photon_energy) {
    double xmin = g_mesh->dx * (double)(node->i - 1),
           xmax = g_mesh->dx * (double)(node->i);

    double alpha = absorption_coefficient(node->material, photon_energy);

    return (int)(g_config->particles_per_cell * 10 * (exp(-alpha * xmin) - exp(-alpha * xmax)));
}


// Calculated relative absorption rates for the different possible transitions at the given
// energy.
int calc_absorption_rates(Material material, double transistion_rate[NOAMTIA][DIME][3]) {
    for(int e = 0; e < DIME; ++e) {
        double energy = (double)e * DE + material.Eg;
        double sum = 0.;

        for(int v = 0; v < 3; ++v) {
            sum += optical_transition_rate(material.id, 1, v, energy);
            transistion_rate[material.id][e][v] = sum;
        }
        for(int v = 0; v < 3; ++v) {
            transistion_rate[material.id][e][v] /= sum;
        }
    }

    return 1;
}


Particle create_photoexcited_carrier(Node *node, double photon_energy, int conduction_band, int valence_band) {
    double mc = node->mat->cb[conduction_band].mstar,
           mv = node->mat->vb[valence_band].mstar;
    double mr = mc * mv / (mc + mv);

    double e2 = (mr / mc) * (photon_energy - node->mat->Eg - node->mat->vb[valence_band].emin);
    printf("%g\n", e2);

    Vec2 loc = mc_random_location_in_node(node);

    double kf = SMH[node->material][conduction_band]
              * sqrt(e2 * (1. + node->mat->cb[conduction_band].alpha * e2));
    double cs  = 1. - 2. * rnd();     // random number -1 -> 1
    double sn  = sqrt(1. - cs * cs);
    double fai = 2. * PI * rnd();     // random angle (radians) 0 -> 2pi
    double kx = kf * cs,
           ky = kf * sn * cos(fai),
           kz = kf * sn * sin(fai);
    int id = g_config->next_particle_id++;
    double time = -log(rnd())/GM[node->material];

    return (Particle){.id=id,
                      .t=time,
                      .valley=conduction_band,
                      .x=loc.x,
                      .y=loc.y,
                      .kx=kx,
                      .ky=ky,
                      .kz=kz};
}


int photoexcite_carriers(Mesh *mesh, double photon_energy, double transistion_rate[NOAMTIA][DIME][3]) {

    int p = 0;
    for(int j = 1; j <= mesh->ny + 1; ++j) {
        for(int i = 1; i <= mesh->nx + 1; ++i) {
            Node *node = mc_node(i, j);
            int e = (int)((photon_energy - node->mat->Eg) / DE);
            int num_carriers = electrons_in_cell(node, photon_energy);

            for(int n = 0; n < num_carriers; ++n) {
                double r = rnd();
                for(int v = 0; v < 3; ++v) {
                    if(r <= transistion_rate[node->material][e][v]) {
                        P[p] = create_photoexcited_carrier(node, photon_energy, 1, v);
                        ++p;
                        break;
                    }
                }
            } // carriers per cell

        } // i
    } // j
    g_config->num_particles = p;

    return 0;
}
