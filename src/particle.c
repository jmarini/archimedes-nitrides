#include "particle.h"

#include "configuration.h"
#include "global_defines.h"
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
    return g_config->next_particle_id++;
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
