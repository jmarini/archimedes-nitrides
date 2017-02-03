#include "particle.h"

#include "configuration.h"
#include "constants.h"
#include "global_defines.h"
#include "random.h"
#include "vec.h"


Index mc_particle_coords(Particle *particle) {
    // uses globabl variables dx & dy
    int i = (int)(particle->x / g_mesh->dx) + 1;
    if(i < 1) { i = 1; }
    if(i > g_mesh->nx ) { i = g_mesh->nx; }

    int j = (int)(particle->y / g_mesh->dy) + 1;
    if(j < 1) { j = 1; }
    if(j > g_mesh->ny ) { j = g_mesh->ny; }

    return (Index){.i=i, .j=j};
}


Index mc_particle_edge_coords(Particle *particle) {
    // uses globabl variables dx & dy
    int i = (int)(particle->x / g_mesh->dx + 1.5);
    int j = (int)(particle->y / g_mesh->dy + 1.5);

    return (Index){.i=i, .j=j};
}


Node * mc_get_particle_node(Particle *particle) {
    return mc_node_s(mc_particle_coords(particle));
}


long long int mc_next_particle_id( ) {
    static long long int next_particle_id = 0;
    return next_particle_id++;
}


double mc_particle_energy(Particle *particle) {
    Node *node = mc_get_particle_node(particle);

    if(g_config->conduction_band == PARABOLIC) {
        return node->mat->cb.hhm[particle->valley] * mc_particle_ksquared(particle);
    }
    else if(g_config->conduction_band == KANE) {
        real alpha = node->mat->cb.alpha[particle->valley];
        real gamma = node->mat->cb.hhm[particle->valley] * mc_particle_ksquared(particle);
        return (-1.0 + sqrt(1.0 + 4.0 * alpha * gamma)) / (2.0 * alpha);
    }
    else {
        return -1.0;
    }
}


double mc_particle_norm_energy(Particle *particle, int axis) {
    Node *node = mc_get_particle_node(particle);

    double ksquared = 0.;
    switch(axis) {
        case 0: // x
            ksquared = particle->kx * particle->kx; break;
        case 1: // y
            ksquared = particle->ky * particle->ky; break;
        case 2: // z
            ksquared = particle->kz * particle->kz; break;
        default:
            ksquared = mc_particle_ksquared(particle);
    }

    if(g_config->conduction_band == PARABOLIC) {
        return node->mat->cb.hhm[particle->valley] * ksquared;
    }
    else if(g_config->conduction_band == KANE) {
        real alpha = node->mat->cb.alpha[particle->valley];
        real gamma = node->mat->cb.hhm[particle->valley] * ksquared;
        return (-1.0 + sqrt(1.0 + 4.0 * alpha * gamma)) / (2.0 * alpha);
    }
    else {
        return -1.0;
    }
}


int mc_calculate_isotropic_k(Particle *particle, double new_energy) {
    double k = 0.;
    Node *node = mc_get_particle_node(particle);
    if(g_config->conduction_band == KANE) {
        k = node->mat->cb.smh[particle->valley]
          * sqrt(new_energy * (1. + node->mat->cb.alpha[particle->valley] * new_energy));
    }
    else if(g_config->conduction_band == PARABOLIC) {
        k = node->mat->cb.smh[particle->valley] * sqrt(new_energy);
    }
    else { return 0; }

    double cs  = 1. - 2. * rnd( );
    double sn  = sqrt(1. - cs * cs);
    double fai = 2. * PI * rnd( );

    particle->kx = k * cs;
    particle->ky = k * sn * cos(fai);
    particle->kz = k * sn * sin(fai);

    return 1;
}
