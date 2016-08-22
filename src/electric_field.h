/* electric_field.h -- This file is part of GNU archimedes

   Archimedes is a simulator for Submicron and Nanoscaled
   2D III-V Semiconductor Devices.

   Copyright (C) 2004-2011 Jean Michel D. Sellier
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


// Computation of the electrostatic potential,
// i.e. resolution of the 2D Poisson equation,
// by means of the computation of the stationary
// solution of a pseudo-transient Poisson equation.
// From version 0.1.0 on, the potential is added to the
// the minimum energy of the semiconductor material
// in order to take into account heterostructures.
// For more information see the manual of
// GNU Archimedes release 1.0.0.


void Electric_Field(void) {
    int i = 0,
        j = 0;
    int nx = g_mesh->nx,
        ny = g_mesh->ny;
    real dx = g_mesh->dx,
         dy = g_mesh->dy;
    real factor = 0.9,
         kappa  = 0.,
         deltat = 0.,
         rho    = 0.;
    real dx2 = 1. / (dx * dx),
         dy2 = 1. / (dy * dy);
    real potential[NXM + 1][NYM + 1];
    real neighbors_x = 0.,
         neighbors_y = 0.;

    if(SIO2_UP_FLAG || SIO2_DOWN_FLAG) {
        printf("Error: SIO2 flag is deprecated.\n");
    }

    poisson_boundary_conditions( );

    for(int n = 0; n < POISSONITMAX; ++n) {
        poisson_boundary_conditions( );

        for(i = 1; i <= nx + 1; ++i) {
            for(j = 1; j <= ny + 1; ++j) {
                potential[i][j] = g_mesh->info[i][j].potential;
            }
        }

        // exclude edge nodes in iteration
        for(j = 2; j <= ny; ++j) {
            for(i = 2; i <= nx; ++i) {
                kappa = EPSR[g_mesh->info[i][j].material] * EPS0 / Q;
                deltat = factor * 0.5 / (kappa * (dx2 + dy2));
                rho = (g_mesh->info[i][j].e.density - g_mesh->info[i][j].donor_conc)
                    - (g_mesh->info[i][j].h.density - g_mesh->info[i][j].acceptor_conc); // charge neutrality eqn.

                // here we are calculating the difference in potential
                // between nearest neighbors
                //   e.g. for x-axis: (p(i+1,j) - p(i,j)) - (p(i,j) - p(i-1,j))
                neighbors_x = potential[i+1][j  ] - 2. * potential[i][j] + potential[i-1][j  ];
                neighbors_y = potential[i  ][j+1] - 2. * potential[i][j] + potential[i  ][j-1];
                g_mesh->info[i][j].potential = potential[i][j]
                                             - deltat * rho
                                             + deltat * kappa * (neighbors_x * dx2 + neighbors_y * dy2);
            }
        }
    }

    poisson_boundary_conditions( );

    // We save the classical potential
    // and we substract the energy minimum of the
    // semiconductor material in order to
    // take into account heterostructures
    for(i = 1; i <= nx + 1; ++i) {
        for(j = 1; j <= ny + 1; ++j) {
            g_mesh->info[i][j].potential -= EMIN[g_mesh->info[i][j].material][1];
        }
    }

    if(g_config->quantum_flag){
        printf("Calculation of Quantum Effective Potential\n");
        // We take in account the Quantum Effects
        quantum_effective_potential( );
        poisson_boundary_conditions( );
    }

    // Computation of the X-component of the electric Field
    // ====================================================
    for(j = 1; j <= ny + 1; ++j) {
        for(i = 2; i <= nx; ++i) { // calculate e-field at edges separately
            g_mesh->info[i][j].efield_x =
                -0.5 * (g_mesh->info[i+1][j].potential - g_mesh->info[i-1][j].potential) / dx;
        }

        // set electric field at edges
        g_mesh->info[     1][j].efield_x = g_mesh->info[ 2][j].efield_x;
        g_mesh->info[nx + 1][j].efield_x = g_mesh->info[nx][j].efield_x;
    }

    // Computation of the Y-component of the electric Field
    // ====================================================
    for(i = 1; i <= nx + 1; ++i) {
        for(j = 2; j <= ny; ++j) {
            g_mesh->info[i][j].efield_y =
                -0.5 * (g_mesh->info[i][j+1].potential - g_mesh->info[i][j-1].potential) / dy;
        }

        // set electric field at edges
        g_mesh->info[i][     1].efield_y = g_mesh->info[i][ 2].efield_y;
        g_mesh->info[i][ny + 1].efield_y = g_mesh->info[i][ny].efield_y;
    }

    // compatability
    for(i = 1; i <= nx + 1; ++i) {
        for(j = 1; j <= ny + 1; ++j) {
            E[i][j][0] = g_mesh->info[i][j].efield_x;
            E[i][j][1] = g_mesh->info[i][j].efield_y;
        }
    }

}
