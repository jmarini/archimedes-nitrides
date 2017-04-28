#include "particle_creation.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "configuration.h"
#include "constants.h"
#include "mesh.h"
#include "particle.h"
#include "random.h"
#include "vec.h"




/* Calculate number of superparticles for the specified node
 */
static int superparticles_per_cell(Mesh *mesh, Node *node) {
    // Load data from previous simulation
    if(g_config->load_initial_data == ON || g_config->tcad_data == ON) {
        return (int)ceil(node->e.density * mesh->dx * mesh->dy
                         / g_config->carriers_per_superparticle);
    }
    // Populate from doping
    else {
        return (int)ceil(node->donor_conc * mesh->dx * mesh->dy
                         / g_config->carriers_per_superparticle);
    }
}


/* Select energy based on either
 */
static double select_energy(Node *node) {
    double r = 0.;
    if(g_config->load_initial_data == ON) {
        r = node->e.energy / node->e.density; // TODO: why??
    }
    else {
        r = -log(rnd( )) * KB * g_config->lattice_temp; // TODO: why negative?
    }
    return 1.5 * r / Q;
}


/* Select a random valley, upper is percent chance to be in satellite valley
 */
static int select_valley(Material *material, double upper) {
    if(material->cb.num_valleys > 1 && rnd( ) > upper) {
        return 2;
    }

    return 1;
}


/* Choose random isotropic k vector for given energy and valley
   1. Calculate |k|
   2. Select random cos(theta)
   3. Select random phi
 */
static Vec3 select_isotropic_k(Material *material, double energy, int valley) {
    double k = 0.;
    if(g_config->conduction_band == KANE) {
        k = material->cb.smh[valley]
          * sqrt(energy * (1. + material->cb.alpha[valley] * energy));
    }
    else if(g_config->conduction_band == PARABOLIC) {
        k = material->cb.smh[valley] * sqrt(energy);
    }
    else { return (Vec3){.x=0., .y=0., .z=0.}; }

    double costheta = 1. - 2. * rnd( );
    double sintheta = sqrt(1. - costheta * costheta);
    double phi = 2. * PI * rnd( );


    return (Vec3){.x=k * costheta * sin(phi),
                  .y=k * sintheta * sin(phi),
                  .z=k * cos(phi)};
}



/* Choose random isotropic k vector for given energy and valley
   1. Calculate |k|
   2. Select random cos(theta)
   3. Select random phi
 */
static Vec3 select_isotropic_k_edge(Material *material, double energy, int valley, int direction) {
    double k = 0.;
    if(g_config->conduction_band == KANE) {
        k = material->cb.smh[valley]
          * sqrt(energy * (1. + material->cb.alpha[valley] * energy));
    }
    else if(g_config->conduction_band == PARABOLIC) {
        k = material->cb.smh[valley] * sqrt(energy);
    }
    else { return (Vec3){.x=0., .y=0., .z=0.}; }

    double costheta = rnd( ); // theta between 0 and pi/2
    double sintheta = sqrt(1. - costheta * costheta);
    double phi = 2. * PI * rnd( );

    Vec2 sign = {1., 1.};

    if(direction == direction_t.TOP) {        sign.y *= -1.; }
    else if(direction == direction_t.RIGHT) { sign.x *= -1.; }

    return (Vec3){.x=sign.x * k * costheta,
                  .y=sign.y * k * sintheta * sin(phi),
                  .z=k * sintheta * cos(phi)};
}


/* Select a random position in the node, offset by dx/2
   If the node is at an edge, random position within dx/2
   If node is in middle, random position in full node, offset by dx/2
 */
static Vec2 select_position(Mesh *mesh, Node *node) {
    double x = mesh->dx * (rnd( ) + (double)(node->i) - 1.5);
    double y = mesh->dy * (rnd( ) + (double)(node->j) - 1.5);

    if(node->i == 1) {
        x = mesh->dx * 0.5 * rnd( );
    }
    if(node->j == 1) {
        y = mesh->dy * 0.5 * rnd( );
    }
    if(node->i == mesh->nx + 1) {
        x = mesh->width - mesh->dx * 0.5 * rnd( );
    }
    if(node->j == mesh->ny + 1) {
        y = mesh->height - mesh->dy * 0.5 * rnd( );
    }

    return (Vec2){.x=x, .y=y};
}


/* Select a random position in the node, offset by dx/2
   If the node is at an edge, random position within dx/2
   If node is in middle, random position in full node, offset by dx/2
 */
static Vec2 select_edge_position(Mesh *mesh, int index, int direction) {
    double x = mesh->dx * (rnd( ) + (double)(index) - 1.5);
    double y = mesh->dy * (rnd( ) + (double)(index) - 1.5);

    if((direction == direction_t.TOP || direction == direction_t.BOTTOM) &&  index == 1) {
        x = mesh->dx * 0.5 * rnd( );
    }
    if((direction == direction_t.TOP || direction == direction_t.BOTTOM) &&  index == mesh->nx) {
        x = mesh->width - mesh->dx * 0.5 * rnd( );
    }
    if(direction == direction_t.BOTTOM) { y = mesh->dy * 0.5 * rnd(); }
    if(direction == direction_t.TOP) { y = mesh->height - mesh->dy * 0.5 * rnd(); }

    if((direction == direction_t.LEFT || direction == direction_t.RIGHT) &&  index == 1) {
        y = mesh->dy * 0.5 * rnd( );
    }
    if((direction == direction_t.LEFT || direction == direction_t.RIGHT) &&  index == mesh->ny) {
        y = mesh->height - mesh->dy * 0.5 * rnd( );
    }
    if(direction == direction_t.LEFT) { x = mesh->dx * 0.5 * rnd(); }
    if(direction == direction_t.RIGHT) { x = mesh->width - mesh->dx * 0.5 * rnd(); }

    return (Vec2){.x=x, .y=y};
}


/* Select random time based on scattering rates
 */
static double select_time(Material *material, double total_scattering_rate[NOAMTIA+1]) {
    return -log(rnd( )) / total_scattering_rate[material->id];
}



/* Create particle at node
 */
Particle create_particle(Mesh *mesh, Node *node, double upper_valley, double total_scattering_rate[NOAMTIA+1]) {
    double energy = select_energy(node);
    int valley = select_valley(node->material, upper_valley);
    Vec3 k = select_isotropic_k(node->material, energy, valley);
    double time = select_time(node->material, total_scattering_rate);
    Vec2 pos = select_position(mesh, node);

    return (Particle){
        .id=mc_next_particle_id( ),
        .valley=valley,
        .kx=k.x,
        .ky=k.y,
        .kz=k.z,
        .x=pos.x,
        .y=pos.y,
        .t=time
    };
}


/* Create particle at edge
 */
Particle create_edge_particle(Mesh *mesh, int index, int direction, double start_time, double upper_valley, double total_scattering_rate[NOAMTIA+1]) {
    Vec2 pos = select_edge_position(mesh, index, direction);
    int i = (int)(pos.x / mesh->dx + 1.5);
    int j = (int)(pos.y / mesh->dy + 1.5);
    if(i <= 1) { i = 1; }
    if(j <= 1) { j = 1; }
    if(i >= mesh->nx + 1) { i = mesh->nx + 1; }
    if(j >= mesh->ny + 1) { j = mesh->ny + 1; }
    Node *node = mc_node(i, j);

    int valley = select_valley(node->material, upper_valley);
    double energy = select_energy(node);
    Vec3 k = select_isotropic_k_edge(node->material, energy, valley, direction);
    double time = start_time + select_time(node->material, total_scattering_rate);

    return (Particle){
        .id=mc_next_particle_id( ),
        .valley=valley,
        .kx=k.x,
        .ky=k.y,
        .kz=k.z,
        .x=pos.x,
        .y=pos.y,
        .t=time
    };
}


/*  For each node, calculate the number of superparticles depending on the
    specified models. Then create that many particles assuming an isotropic
    distribution of k.
 */
int populate_superparticles(Mesh *mesh, double upper_valley, double total_scattering_rate[NOAMTIA+1]) {
    int index = 0;

    for(int i = 1; i <= mesh->nx + 1; ++i) {
        for(int j = 1; j <= mesh->ny + 1; ++j) {
            Node *node = &(mesh->nodes[i][j]);

            int sppc = superparticles_per_cell(mesh, node);
            if((i == 1) || (i == mesh->nx + 1)) { sppc /= 2; }
            if((j == 1) || (j == mesh->ny + 1)) { sppc /= 2; }

            if(index + sppc > NPMAX) {
                printf("ERROR: Number of particles exceeds maximum (%d)\n", NPMAX);
                exit(EXIT_FAILURE);
            }

            for(int n = 0; n < sppc; ++n) {
                ++index;
                mesh->particles[index] = create_particle(mesh, node, upper_valley, total_scattering_rate);
            }
        }
    }

    g_config->num_particles = index;

    return 0;
}
