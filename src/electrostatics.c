#include "electrostatics.h"

#include <math.h>
#include <stdio.h>

#include "configuration.h"
#include "constants.h"
#include "global_defines.h"
#include "mesh.h"


// =============================
// === Convenience Functions ===
// =============================

// Convenience function to calculate potential and electric/magnetic fields
int electrostatics(Mesh *mesh) {
    int error = 0;

    if(g_config->poisson_flag == ON) { error |= poisson(mesh); }
    if(g_config->faraday_flag == ON) { error |= faraday(mesh); }

    if(error != 0) {
        printf("Error: Unknown error solving electrostatics.\n");
        return 1;
    }

    return 0;
}


// Convenience function to calculate potential and electric field
int poisson(Mesh *mesh) {
    // Classical potential
    if(calculate_potential(mesh) != 0) {
        printf("Error: Unknown error calculating potential.\n");
        return 1;
    }

    // Electric field
    if(electric_field(mesh) != 0) {
        printf("Error: Unknown error calculating electric field.\n");
        return 1;
    }

    return 0;
}


// Convenience function to calculate boundary conditions and magnetic field
int faraday(Mesh *mesh) {
    // Faraday boundary conditions
    if(faraday_boundary_conditions(mesh) != 0) {
        printf("Error: Unknown error calculating Faraday boundary conditions.\n");
        return 1;
    }

    // Mangetic Field
    if(magnetic_field(mesh) != 0) {
        printf("Error: Unknown error calculating magnetic field.\n");
        return 1;
    }

    return 0;
}


// ===========================
// === Boundary Conditions ===
// ===========================

int faraday_boundary_conditions(Mesh *mesh) {
    int nx = mesh->nx,
        ny = mesh->ny;

    for(int i = 1; i <= nx + 1; ++i) {
        // Bottom Edge
        mesh->nodes[i][0].magnetic_field = mesh->nodes[i][3].magnetic_field;
        mesh->nodes[i][1].magnetic_field = mesh->nodes[i][2].magnetic_field;

        // Upper Edge
        mesh->nodes[i][ny+1].magnetic_field = mesh->nodes[i][ny  ].magnetic_field;
        mesh->nodes[i][ny+2].magnetic_field = mesh->nodes[i][ny-1].magnetic_field;
    }

    for(int j = 1; j <= ny + 1; ++j) {
        // Left Edge
        mesh->nodes[0][j].magnetic_field = mesh->nodes[3][j].magnetic_field;
        mesh->nodes[1][j].magnetic_field = mesh->nodes[2][j].magnetic_field;

        // Right Edge
        mesh->nodes[nx+1][j].magnetic_field = mesh->nodes[nx-1][j].magnetic_field;
        mesh->nodes[nx+2][j].magnetic_field = mesh->nodes[nx  ][j].magnetic_field;
    }

    return 0;
}


int poisson_boundary_conditions(Mesh *mesh) {
    int nx = mesh->nx,
        ny = mesh->ny;

    // These are completely generic boundary conditions
    for(int i = 1; i <= nx + 1; ++i) {
        // Bottom Edge
        // ===========

        // insulator
        if(mc_is_boundary_insulator(direction_t.BOTTOM, i) || mc_is_boundary_vacuum(direction_t.BOTTOM, i)) {
            // if there is no potential on insulator
            mesh->nodes[i][1].potential = mesh->nodes[i][2].potential;

            // if there is potential on the insulator
            //   set the potential on the first two nodes
            if(mesh->edges[direction_t.BOTTOM][i].potential != 0.0) {
                mesh->nodes[i][1].potential = mesh->edges[direction_t.BOTTOM][i].potential;
            }
        }
        // schottky or ohmic contact
        //   set potential of first two nodes to that of contact
        else if(mc_is_boundary_contact(direction_t.BOTTOM, i)) {
            mesh->nodes[i][1].potential = mesh->edges[direction_t.BOTTOM][i].potential;
        }

        // Upper Edge
        // ==========

        // insulator
        if(mc_is_boundary_insulator(direction_t.TOP, i) || mc_is_boundary_vacuum(direction_t.TOP, i)) {
            // if there is no potential on insulator
            mesh->nodes[i][ny + 1].potential = mesh->nodes[i][ny    ].potential;

            // if there is potential on the insulator
            //   set the potential on the first two nodes
            if(mesh->edges[direction_t.TOP][i].potential != 0.0) {
                mesh->nodes[i][ny + 1].potential = mesh->edges[direction_t.TOP][i].potential;
            }
        }
        // schottky or ohmic contact
        //   set potential of first two nodes to that of contact
        else if(mc_is_boundary_contact(direction_t.TOP, i)) {
            mesh->nodes[i][ny + 1].potential = mesh->edges[direction_t.TOP][i].potential;
        }

    }


    for(int j = 1; j <= ny + 1; ++j) {
        // Left Edge
        // =========

        // insulator
        if(mc_is_boundary_insulator(direction_t.LEFT, j) || mc_is_boundary_vacuum(direction_t.LEFT, j)) {
            mesh->nodes[1][j].potential = mesh->nodes[2][j].potential;

            // if there is potential on the insulator
            //   set the potential on the first two nodes
            if(mesh->edges[direction_t.LEFT][j].potential != 0.0) {
                mesh->nodes[1][j].potential = mesh->edges[direction_t.LEFT][j].potential;
            }
        }
        // schottky or ohmic contact
        //   set potential of first two nodes to that of contact
        else if(mc_is_boundary_contact(direction_t.LEFT, j)) {
            mesh->nodes[1][j].potential = mesh->edges[direction_t.LEFT][j].potential;
        }

        // Right Edge
        // ==========

        // insulator
        if(mc_is_boundary_insulator(direction_t.RIGHT, j) || mc_is_boundary_vacuum(direction_t.RIGHT, j)) {
            mesh->nodes[nx + 1][j].potential = mesh->nodes[nx - 1][j].potential;

            // if there is potential on the insulator
            //   set the potential on the first two nodes
            if(mesh->edges[direction_t.RIGHT][j].potential != 0.0) {
                mesh->nodes[nx + 1][j].potential = mesh->edges[direction_t.RIGHT][j].potential;
            }
        }
        // schottky or ohmic contact
        //   set potential of first two nodes to that of contact
        else if(mc_is_boundary_contact(direction_t.RIGHT, j)) {
            mesh->nodes[nx + 1][j].potential = mesh->edges[direction_t.RIGHT][j].potential;
        }

    }

    return 0;

}


// ==========================
// === Potential / Fields ===
// ==========================


int calculate_potential(Mesh *mesh) {
    int nx = mesh->nx,
        ny = mesh->ny;

    real factor = 0.9; // successive over-relaxation factor

    if(poisson_boundary_conditions(mesh) != 0) {
        printf("Error: Unknown error calculating Poisson boundary conditions.\n");
        return 1;
    }

    double potential[NXM + 1][NYM + 1]; // holds potential from previous iteration
    // TODO: stop iteration based on residual OR iteration max - should lead to faster sim
    for(int n = 0; n < POISSONITMAX; ++n) {
        if(poisson_boundary_conditions(mesh) != 0) {
            printf("Error: Unknown error calculating Poisson boundary conditions.\n");
            return 1;
        }

        // potential from previous iteration
        for(int i = 1; i <= nx + 1; ++i) {
            for(int j = 1; j <= ny + 1; ++j) {
                potential[i][j] = mesh->nodes[i][j].potential;
            }
        }

        // exclude edge nodes
        for(int j = 2; j <= ny; ++j) {
            for(int i = 2; i <= nx; ++i) {
                double dx2 = mesh->dx * mesh->dx,
                       dy2 = mesh->dy * mesh->dy;
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

    if(poisson_boundary_conditions(mesh) != 0) {
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
        quantum_effective_potential(mesh);
        if(poisson_boundary_conditions(mesh) != 0) {
            printf("Error: Unknown error calculating Poisson boundary conditions.\n");
            return 1;
        }
    }

    return 0;
}


int quantum_effective_potential(Mesh *mesh) {
    return 1; // not implemented
}


// Calculate the electric field for a given calculated potential
int electric_field(Mesh *mesh) {
    int nx = mesh->nx,
        ny = mesh->ny;

    // X-component of the electric field
    for(int j = 1; j <= ny + 1; ++j) {
        for(int i = 2; i <= nx; ++i) { // calculate e-field at edges separately
            double deltaV = mesh->nodes[i+1][j].potential - mesh->nodes[i-1][j].potential;
            mesh->nodes[i][j].efield.x = -0.5 * deltaV / mesh->dx;
        }

        // set electric field at edges
        mesh->nodes[     1][j].efield.x = mesh->nodes[ 2][j].efield.x;
        mesh->nodes[nx + 1][j].efield.x = mesh->nodes[nx][j].efield.x;
    }

    // Y-component of the electric Field
    for(int i = 1; i <= nx + 1; ++i) {
        for(int j = 2; j <= ny; ++j) {
            double deltaV = mesh->nodes[i][j+1].potential - mesh->nodes[i][j-1].potential;
            mesh->nodes[i][j].efield.y = -0.5 * deltaV / mesh->dy;
        }

        // set electric field at edges
        mesh->nodes[i][     1].efield.y = mesh->nodes[i][ 2].efield.y;
        mesh->nodes[i][ny + 1].efield.y = mesh->nodes[i][ny].efield.y;
    }

    return 0;
}


int magnetic_field(Mesh *mesh) {
    for(int i = 2; i <= mesh->nx; ++i) {
        for(int j = 2; j <= mesh->ny; ++j) {
            double delEx = 0.5 * (mesh->nodes[i  ][j+1].efield.x - mesh->nodes[i  ][j-1].efield.x) / mesh->dy;
            double delEy = 0.5 * (mesh->nodes[i+1][  j].efield.y - mesh->nodes[i-1][  j].efield.y) / mesh->dx;

            mesh->nodes[i][j].magnetic_field = 0.25
                * (mesh->nodes[i+1][j].magnetic_field + mesh->nodes[i][j+1].magnetic_field +
                   mesh->nodes[i-1][j].magnetic_field + mesh->nodes[i][j-1].magnetic_field)
                - g_config->dt * (delEy - delEx);
        }
    }

    return 0;
}
