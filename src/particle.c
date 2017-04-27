#include "particle.h"

#include <stdio.h>

#include "configuration.h"
#include "constants.h"
#include "global_defines.h"
#include "material.h"
#include "random.h"
#include "vec.h"



static inline int clamp(int value, int min, int max) {
    int ret = value < min ? min : value;
    return ret > max ? max : ret;
}


Index mc_particle_coords(Particle *p) {
    int i = clamp((int)(p->x / g_mesh->dx) + 1, 1, g_mesh->nx);
    int j = clamp((int)(p->y / g_mesh->dy) + 1, 1, g_mesh->ny);

    return (Index){.i=i, .j=j};
}


Index mc_particle_edge_coords(Particle *p) {
    int i = (int)(p->x / g_mesh->dx + 1.5);
    int j = (int)(p->y / g_mesh->dy + 1.5);

    return (Index){.i=i, .j=j};
}


Node * mc_get_particle_node(Particle *p) {
    return mc_node_s(mc_particle_coords(p));
}


long long int mc_next_particle_id( ) {
    static long long int next_particle_id = 0;
    return next_particle_id++;
}


double mc_particle_energy(Particle *p) {
    Material *material = mc_get_particle_node(p)->material;

    if(g_config->conduction_band == PARABOLIC) {
        return material->cb.hhm[p->valley] * mc_particle_ksquared(p);
    }
    else if(g_config->conduction_band == KANE) {
        real alpha = material->cb.alpha[p->valley];
        real gamma = material->cb.hhm[p->valley] * mc_particle_ksquared(p);
        return (-1.0 + sqrt(1.0 + 4.0 * alpha * gamma)) / (2.0 * alpha);
    }
    else {
        return -1.0;
    }
}


double mc_particle_ksquared(Particle *p) {
    return (p->kx * p->kx +
            p->ky * p->ky +
            p->kz * p->kz);
}


double mc_particle_k(Particle *p) {
    return sqrt(mc_particle_ksquared(p));
}


double mc_particle_norm_energy(Particle *p, int axis) {
    Material *material = mc_get_particle_node(p)->material;

    double ksquared = 0.;
    switch(axis) {
        case 0: // x
            ksquared = p->kx * p->kx; break;
        case 1: // y
            ksquared = p->ky * p->ky; break;
        case 2: // z
            ksquared = p->kz * p->kz; break;
        default:
            ksquared = mc_particle_ksquared(p);
    }

    if(g_config->conduction_band == PARABOLIC) {
        return material->cb.hhm[p->valley] * ksquared;
    }
    else if(g_config->conduction_band == KANE) {
        real alpha = material->cb.alpha[p->valley];
        real gamma = material->cb.hhm[p->valley] * ksquared;
        return (-1.0 + sqrt(1.0 + 4.0 * alpha * gamma)) / (2.0 * alpha);
    }
    else {
        return -1.0;
    }
}


int mc_calculate_isotropic_k(Particle *p, double new_energy) {
    Material *material = mc_get_particle_node(p)->material;

    double k = 0.;
    if(g_config->conduction_band == KANE) {
        k = material->cb.smh[p->valley]
          * sqrt(new_energy * (1. + material->cb.alpha[p->valley] * new_energy));
    }
    else if(g_config->conduction_band == PARABOLIC) {
        k = material->cb.smh[p->valley] * sqrt(new_energy);
    }
    else { return 1; }

    double cs  = 1. - 2. * rnd( );
    double sn  = sqrt(1. - cs * cs);
    double fai = 2. * PI * rnd( );

    p->kx = k * cs;
    p->ky = k * sn * cos(fai);
    p->kz = k * sn * sin(fai);

    return 0;
}


int mc_calculate_anisotropic_k(Particle *p, double ki, double kf, double cb) {
    double sb  = sqrt(1. - cb * cb);
    double fai = 2. * PI * rnd();
    double skk = sqrt(p->kx * p->kx +
                      p->ky * p->ky);
    double a11 =  p->ky / skk;
    double a12 =  p->kx * p->kz / skk / ki;
    double a13 =  p->kx / ki;
    double a21 = -p->kx / skk;
    double a22 =  p->ky * p->kz / skk / ki;
    double a23 =  p->ky / ki;
    double a32 = -skk / ki;
    double a33 =  p->kz / ki;
    double x1 = kf * sb * cos(fai);
    double x2 = kf * sb * sin(fai);
    double x3 = kf * cb;

    p->kx = a11 * x1 + a12 * x2 + a13 * x3;
    p->ky = a21 * x1 + a22 * x2 + a23 * x3;
    p->kz =            a32 * x2 + a33 * x3;

    return 0;
}


particle_info_t mc_calculate_particle_info(Particle *p) {
    // calculate particle coordinates
    int i = 0,
        j = 0;
    int nx = g_mesh->nx,
        ny = g_mesh->ny;
    double dx = g_mesh->dx,
           dy = g_mesh->dy;

    i = (int)(p->x / dx + 1.5);
    if(i <= 1) { i = 1; }
    if(i >= nx + 1) { i = nx + 1; }

    j = (int)(p->y / dy + 1.5);
    if(j <= 1) { j = 1; }
    if(j >= ny + 1) { j = ny + 1; }

    Material *material = g_mesh->nodes[i][j].material;

    // calculate particle energy and velocity
    double ksquared = mc_particle_ksquared(p);
    double energy = 0.;
    double xvelocity = 0.,
           yvelocity = 0.;

    if(g_config->conduction_band == PARABOLIC) {
        energy = material->cb.hhm[p->valley] * ksquared;
        xvelocity = p->kx * material->cb.hm[p->valley];
        yvelocity = p->ky * material->cb.hm[p->valley];
    }
    else if(g_config->conduction_band == KANE) {
        double sq = sqrt(1. + 4. * material->cb.alpha[p->valley]
                                 * material->cb.hhm[p->valley] * ksquared);
        energy = (sq - 1.) / (2. * material->cb.alpha[p->valley]);
        xvelocity = p->kx * material->cb.hm[p->valley] / sq;
        yvelocity = p->ky * material->cb.hm[p->valley] / sq;
    }


    return (particle_info_t){
        .id=p->id,
        .valley=p->valley,
        .kx=p->kx,
        .ky=p->ky,
        .kz=p->kz,
        .energy=energy,
        .t=p->t,
        .x=p->x,
        .y=p->y,
        .i=i,
        .j=j,
        .vx=xvelocity,
        .vy=yvelocity
    };
}


void mc_print_tracking(int it, Particle *p) {
    particle_tracking_t pt = (particle_tracking_t){
        .time=(float)p->t,
        .x=(float)p->x,
        .y=(float)p->y,
        .energy=(float)mc_particle_energy(p),
        .valley=(int)p->valley
    };

    char s[150];
    sprintf(s, "tracking%06lld.bin", p->id);

    FILE *fp = NULL;
    if(it == 0.) { fp = fopen(s, "wb"); }
    else {         fp = fopen(s, "ab"); }

    fwrite(&pt, sizeof(particle_tracking_t), 1, fp);
    fclose(fp);
}
