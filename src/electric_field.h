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


int electric_field(void) {
    int nx = g_mesh->nx,
        ny = g_mesh->ny;
    real dx = g_mesh->dx,
         dy = g_mesh->dy;
    real factor = 0.9; // successive over-relaxation factor

    if(poisson_boundary_conditions( ) != 0) {
        printf("Error: Unknown error calculating Poisson boundary conditions.\n");
        return 1;
    }

    real potential[NXM + 1][NYM + 1];
    // TODO: stop iteration based on residual OR iteration max - should lead to faster sim
    for(int n = 0; n < POISSONITMAX; ++n) {
        if(poisson_boundary_conditions( ) != 0) {
            printf("Error: Unknown error calculating Poisson boundary conditions.\n");
            return 1;
        }

        // potential from previous iteration
        for(int i = 1; i <= nx + 1; ++i) {
            for(int j = 1; j <= ny + 1; ++j) {
                potential[i][j] = g_mesh->nodes[i][j].potential;
            }
        }

        // exclude edge nodes
        for(int j = 2; j <= ny; ++j) {
            for(int i = 2; i <= nx; ++i) {
                real dx2 = dx * dx,
                     dy2 = dy * dy;
                Node *node = mc_node(i, j);
                real kappa = node->mat->eps_static * EPS0 / Q;
                real deltat = (factor / kappa) * (dx2 * dy2)
                            / ((2 * (dx2 + dy2)) + dx2 * dy2);
                real rho = (node->e.density - node->donor_conc)
                         - (node->h.density - node->acceptor_conc); // charge neutrality eqn.

                // here we are calculating the difference in potential
                // between nearest neighbors
                //   e.g. for x-axis: (p(i+1,j) - p(i,j)) - (p(i,j) - p(i-1,j))
                real neighbors_x = potential[i+1][j  ] - 2. * potential[i][j] + potential[i-1][j  ];
                real neighbors_y = potential[i  ][j+1] - 2. * potential[i][j] + potential[i  ][j-1];
                node->potential = potential[i][j]
                                - deltat * rho
                                + deltat * kappa * (neighbors_x / dx2 + neighbors_y / dy2);
            }
        }
    }

    if(poisson_boundary_conditions( ) != 0) {
        printf("Error: Unknown error calculating Poisson boundary conditions.\n");
        return 1;
    }
    // We save the classical potential and we subtract the energy minimum of the
    // semiconductor material in order to take into account heterostructures
    for(int i = 1; i <= nx + 1; ++i) {
        for(int j = 1; j <= ny + 1; ++j) {
            Node *node = mc_node(i, j);
            node->potential -= node->mat->cb.emin[1];
        }
    }

    if(g_config->quantum_flag){
        printf("Calculation of Quantum Effective Potential\n");
        // We take in account the Quantum Effects
        quantum_effective_potential( );
        if(poisson_boundary_conditions( ) != 0) {
            printf("Error: Unknown error calculating Poisson boundary conditions.\n");
            return 1;
        }
    }

    // Computation of the X-component of the electric Field
    // ====================================================
    for(int j = 1; j <= ny + 1; ++j) {
        for(int i = 2; i <= nx; ++i) { // calculate e-field at edges separately
            g_mesh->nodes[i][j].efield.x =
                -0.5 * (g_mesh->nodes[i+1][j].potential - g_mesh->nodes[i-1][j].potential) / dx;
        }

        // set electric field at edges
        g_mesh->nodes[     1][j].efield.x = g_mesh->nodes[ 2][j].efield.x;
        g_mesh->nodes[nx + 1][j].efield.x = g_mesh->nodes[nx][j].efield.x;
    }

    // Computation of the Y-component of the electric Field
    // ====================================================
    for(int i = 1; i <= nx + 1; ++i) {
        for(int j = 2; j <= ny; ++j) {
            g_mesh->nodes[i][j].efield.y =
                -0.5 * (g_mesh->nodes[i][j+1].potential - g_mesh->nodes[i][j-1].potential) / dy;
        }

        // set electric field at edges
        g_mesh->nodes[i][     1].efield.y = g_mesh->nodes[i][ 2].efield.y;
        g_mesh->nodes[i][ny + 1].efield.y = g_mesh->nodes[i][ny].efield.y;
    }

    return 0;
}
