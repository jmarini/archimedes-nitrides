#include "optical_absorption.h"

#include "configuration.h"
#include "constants.h"
#include "global_defines.h"
#include "material.h"
#include "mesh.h"
#include "particle.h"
#include "random.h"


// Returns the optical joint density of states for the given conduction and valence bands
// evaluated at the given energy
double optical_joint_DOS(Material material, int conduction_band, int valence_band, double energy) {
    double mc = material.cb.mstar[conduction_band],
           mv = material.vb.mstar[valence_band];
    double mr = mc * mv / (mc + mv);
    double alpha = material.cb.alpha[conduction_band];

    double deltaE = material.cb.emin[conduction_band] + material.Eg;
    double e2 = (mr / mc) * (energy - deltaE);

    double gamma = Q * e2 * (1. + alpha * e2);
    double gamma2 = (1. + 2. * alpha * e2);

    double dos = mr * M * sqrt(2. * mr * M) * sqrt(gamma) * gamma2 / (PI * PI * HBAR * HBAR * HBAR);
    return dos;
}


// Returns the transition rate (1/s) between the given conduction and valence bands at the given
// energy
double optical_transition_rate(Material material, int conduction_band,
                               int valence_band, double photon_energy) {
    if(photon_energy < material.Eg + material.cb.emin[conduction_band] + material.vb.emin[valence_band]) {
        return 0.; // photon energy too small
    }

    double prefactor = PI * HBAR * Q * Q * material.Eg
                     / (3. * material.eps_static * EPS0 * M * photon_energy);
    double rate = prefactor * optical_joint_DOS(material, conduction_band, valence_band, photon_energy);
    return rate;
}


// Returns the absorption coefficient at the given energy. Takes into account direct interband
// transitions only from HH, LH, SO valence bands to G-1 conduction band
double absorption_coefficient(Material material, double photon_energy) {
    double rate = 0.;
    for(int v = 0; v < 3; ++v) {
        rate += optical_transition_rate(material, 1, v, photon_energy) * sqrt(material.eps_static) / VLIGHT;
    }

    return rate * material.abs_correction;
}


// Calculates the number of photoexcited electrons in the given node at the specified energy.
// Photons are assumed to enter device at x=0
// TODO: number of carriers needs to depend on energy - relative W @ E vs max W
int electrons_in_cell(Mesh *mesh, Node *node, double photon_energy) {
    double xmin = mesh->dx * (double)(node->i - 1),
           xmax = mesh->dx * (double)(node->i);

    double alpha = absorption_coefficient(*(node->mat), photon_energy);
    double n = g_config->particles_per_cell / (double)g_mesh->ny;

    return (int)(n * (exp(-alpha * xmin) - exp(-alpha * xmax)));
}


// Calculated relative absorption rates for the different possible transitions at the given
// energy.
int calc_absorption_rates(Material material, double transistion_rate[NOAMTIA][DIME][3]) {
    for(int e = 0; e < DIME; ++e) {
        double energy = (double)e * DE + material.Eg;
        double sum = 0.;

        for(int v = 0; v < 3; ++v) {
            sum += optical_transition_rate(material, 1, v, energy);
            transistion_rate[material.id][e][v] = sum;
        }
        for(int v = 0; v < 3; ++v) {
            transistion_rate[material.id][e][v] /= sum;
        }
    }

    return 1;
}


// Create a photoexcited carrier at the given Node for a the given photon energies.
// The transistion will be between the given cb and vb. Returns the created particle.
// Takes statistical weight as approx the total number of photoexcited carriers
Particle create_photoexcited_carrier(Node *node, double photon_energy,
                                     double total_scattering_rate[NOAMTIA+1],
                                     int conduction_band, int valence_band) {
    Material *material = node->mat;

    double mc = material->cb.mstar[conduction_band],
           mv = material->vb.mstar[valence_band];
    double mr = mc * mv / (mc + mv);

    double e2 = (mr / mc) * (photon_energy - material->Eg - material->vb.emin[valence_band]);

    Vec2 loc = mc_random_location_in_node(node);

    double kf = material->cb.smh[conduction_band]
              * sqrt(e2 * (1. + material->cb.alpha[conduction_band] * e2));
    double cs  = 1. - 2. * rnd();     // random number -1 -> 1
    double sn  = sqrt(1. - cs * cs);
    double fai = 2. * PI * rnd();     // random angle (radians) 0 -> 2pi
    double kx = kf * cs,
           ky = kf * sn * cos(fai),
           kz = kf * sn * sin(fai);
    long long int id = mc_next_particle_id( );
    double time = -log(rnd()) / total_scattering_rate[material->id];

    return (Particle){.id=id,
                      .t=time,
                      .valley=conduction_band,
                      .x=loc.x,
                      .y=loc.y,
                      .kx=kx,
                      .ky=ky,
                      .kz=kz};
}


// Photoexcite carriers at the given energy. Returns the number of photogenerated
// carriers
int photoexcite_carriers(Mesh *mesh, double photon_energy,
                         double transistion_rate[NOAMTIA][DIME][3],
                         double total_scattering_rate[NOAMTIA+1],
                         Particle particles[NPMAX+1]) {
    int p = g_config->num_particles;

    for(int j = 1; j <= mesh->ny + 1; ++j) {
        for(int i = 1; i <= mesh->nx + 1; ++i) {
            Node *node = mc_node(i, j);
            int e = (int)((photon_energy - node->mat->Eg) / DE);
            int num_carriers = electrons_in_cell(mesh, node, photon_energy);

            for(int n = 0; n < num_carriers; ++n) {
                double r = rnd();
                for(int v = 0; v < 3; ++v) {
                    if(r <= transistion_rate[node->material][e][v]) {
                        particles[p] = create_photoexcited_carrier(node, photon_energy, total_scattering_rate, 1, v);
                        if(g_config->tracking_output == ON
                           && particles[p].id % g_config->tracking_mod == 0) {
                          mc_print_tracking(0, &particles[p]);
                        }
                        ++p;
                        break;
                    }
                }
            } // carriers per cell

        } // i
    } // j

    int num_photoexcited = p - g_config->num_particles;
    g_config->num_particles = p;
    return num_photoexcited;
}
