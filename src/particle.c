#include "particle.h"
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
