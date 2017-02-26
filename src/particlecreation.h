/* particlecreation.h -- This file is part of Archimedes release 1.2.0
   Archimedes is a simulator for Submicron 2D III-V semiconductor
   Devices. It implements the Monte Carlo method
   for the simulation of the semiclassical Boltzmann equation for both
   electrons and holes. It also includes the quantum effects by means
   of effective potential method. It is now able to simulate applied
   magnetic fields along with self consistent Faraday equation.

   Copyright (C) 2004, 2005, 2006, 2007 Jean Michel Sellier <sellier@dmi.unict.it>

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
// Created on 06 sep.2004, Siracusa, Italy, J.M.Sellier
// Last modif. : 26 Aug.2011, Carry le Rouet, France, J.M.Sellier
// ######################################################

// create a new super-particle on the edge
// edge = 0 Bottom edge
// edge = 1 Right edge
// edge = 2 Upper edge
// edge = 3 Left edge

inline Particle creation(int i, real t, int edge) {
    int iaux = 0,
        iv   = 0,
        ii   = 0,
        j    = 0;
    real ts = 0.0;
    real kx = 0.0,
         ky = 0.0,
         kz = 0.0,
         x  = 0.0,
         y  = 0.0;
    real c1, c2, c3, c4, c5, c6, c7;
    int nx = g_mesh->nx,
        ny = g_mesh->ny;
    real dx = g_mesh->dx,
         dy = g_mesh->dy;

    // We assume that the particles are initially
    // at near thermal equilibrium

    // creation of the particle position vector r=(X,Y)
    x = dx * (rnd() + (real)i - 1.5);
    y = dy * (rnd() + (real)i - 1.5);
    if((edge == direction_t.BOTTOM || edge == direction_t.TOP) && i == 1) {
        x = dx * 0.5 * rnd();
    }
    if((edge == direction_t.BOTTOM || edge == direction_t.TOP) && i == nx) {
        x = g_mesh->width - dx * 0.5 * rnd();
    }
    if(edge == direction_t.BOTTOM) { y =                  dy * 0.5 * rnd(); }
    if(edge == direction_t.TOP)    { y = g_mesh->height - dy * 0.5 * rnd(); }
    if((edge == direction_t.RIGHT || edge == direction_t.LEFT) && i == 1) {
        y = dy * 0.5 * rnd();
    }
    if((edge == direction_t.RIGHT || edge == direction_t.LEFT) && i == ny) {
        y = g_mesh->height - dy * 0.5 * rnd();
    }
    if(edge == direction_t.RIGHT) { x = g_mesh->width - dx * 0.5 * rnd(); }
    if(edge == direction_t.LEFT)  { x =                 dx * 0.5 * rnd(); }

    // creation of the particle pseudo-wave vector k=(KX,KY,KZ)
    // in the (i,j)-th cell
    ii = (int)(x / dx + 1.5);
    j  = (int)(y / dy + 1.5);
    if(ii <= 1) { ii = 1; }
    if( j <= 1) {  j = 1; }
    if(ii >= nx + 1) { ii = nx + 1; }
    if( j >= ny + 1) {  j = ny + 1; }
    int material = g_mesh->nodes[ii][j].material;
    if(g_materials[material].cb.num_valleys == 1) {
        iv = 1;
        iaux = 0;
    }
    else if(g_materials[material].cb.num_valleys >= 2) {
        iv = iaux = 1;
        // 20% of the created particles belongs to the L-valley
        if(rnd() >= 0.8) { iv = iaux = 2; }
    }
    c1 = log(rnd());
    c2 = g_materials[material].cb.smh[iaux]
       * sqrt(-1.5 * BKTQ * c1 * (1. - g_materials[material].cb.alpha[iv] * 1.5 * BKTQ * c1));
    c3 = rnd();
    c4 = sqrt(1. - c3 * c3);
    c5 = 2. * PI * rnd();
    c6 = sin(c5);
    c7 = cos(c5);
    kx = c2 * c3;
    ky = c2 * c4 * c6;
    kz = c2 * c4 * c7;
    ts = t - log(rnd()) / GM[material];
    if(edge == direction_t.TOP)   { ky *= -1.; }
    if(edge == direction_t.RIGHT) { kx *= -1.; }

    long long int id = mc_next_particle_id( );

    return (Particle){.id=id, .valley=iv, .t=ts, .kx=kx, .ky=ky, .kz=kz, .x=x, .y=y};
}

// =================================================
