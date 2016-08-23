#include "particle.h"
#include "vec.h"


index_s mc_particle_coords(particle_t *particle) {
    // uses globabl variables dx & dy
    int i = (int)(particle->x / g_mesh->dx) + 1;
    if(i < 1) { i = 1; }
    if(i > g_mesh->nx ) { i = g_mesh->nx; }

    int j = (int)(particle->y / g_mesh->dy) + 1;
    if(j < 1) { j = 1; }
    if(j > g_mesh->ny ) { j = g_mesh->ny; }

    return (index_s){.i=i, .j=j};
}


mc_node_t * mc_get_particle_node(particle_t *particle) {
    return mc_node_s(mc_particle_coords(particle));
}
