/* updating.h -- This file is part of GNU archimedes.

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


// Updates the macroscopic variables
// Here we call all the subroutines we
// need for simulating the dynamics of the
// super-particles in the device.
// The solution is writted in the array named u2d.


int updating(int iteration, int model) {
    if(model == MEPE || model == MEPH || model == MEPEH) {
        printf("Error: MEP simulation is deprecated.\n");
        return 1;
    }

    // Computation of the electric field
    // =================================
    if(g_config->poisson_flag == ON) {
        electric_field( );
        if(g_config->faraday_flag) { faraday( ); }
    }


    // Monte Carlo Simulation
    // ======================
    EMC(iteration);
    calculate_particles_per_cell( );
    media( );
    // If timestep would put simulation time after ending time, adjust step
    if(g_config->time + g_config->dt >= g_config->tf) {
        g_config->dt = g_config->tf - g_config->time;
    }
    g_config->time += g_config->dt;


    // Output on some usefull informations about the simulation
    printf("%5d   TIME = %10.4g  (picosec)\n", iteration, g_config->time * 1.e12);
    if(g_config->max_min_output){
        // Compute the maximum and minimum of various macroscopic variables
        int i = 0,
            j = 0;
        real maxi = 0.,
             mini = 0.;

        // Max and Min of Potential
        for(i = 1; i <= g_mesh->nx + 1; ++i) {
            for(j = 1; j <= g_mesh->ny + 1; ++j) {
                if(g_mesh->nodes[i][j].potential >= maxi) { maxi = g_mesh->nodes[i][j].potential; }
                if(g_mesh->nodes[i][j].potential <= mini) { mini = g_mesh->nodes[i][j].potential; }
            }
        }
        printf("Max. Potential = %g V\n", maxi);
        printf("Min. Potential = %g V\n", mini);

        // Max and Min of x-component of electric field
        maxi = mini = 0.;
        for(i = 1; i <= g_mesh-> nx + 1; ++i) {
            for(j = 1; j <= g_mesh-> ny + 1; ++j) {
                if(g_mesh->nodes[i][j].efield.x >= maxi) { maxi = g_mesh->nodes[i][j].efield.x; }
                if(g_mesh->nodes[i][j].efield.x <= mini) { mini = g_mesh->nodes[i][j].efield.x; }
            }
        }
        printf("Max. x-elec.field = %g V/m\n", maxi);
        printf("Min. x-elec.field = %g V/m\n", mini);

        // Max and Min of y-component of electric field
        maxi = mini = 0.;
        for(i = 1; i <= g_mesh-> nx + 1; ++i) {
            for(j = 1; j <= g_mesh-> ny + 1; ++j) {
                if(g_mesh->nodes[i][j].efield.y >= maxi) { maxi = g_mesh->nodes[i][j].efield.y; }
                if(g_mesh->nodes[i][j].efield.y <= mini) { mini = g_mesh->nodes[i][j].efield.y; }
            }
        }
        printf("Max. y-elec.field = %g V/m\n", maxi);
        printf("Min. y-elec.field = %g V/m\n", mini);

        // Max and Min of Density
        maxi = 0.;
        mini = g_config->max_doping;
        for(i = 1; i <= g_mesh->nx + 1; ++i) {
            for(j = 1; j <= g_mesh->ny + 1; ++j) {
                if(g_mesh->nodes[i][j].e.density >= maxi) { maxi = g_mesh->nodes[i][j].e.density; }
                if(g_mesh->nodes[i][j].e.density <= mini) { mini = g_mesh->nodes[i][j].e.density; }
            }
        }
        printf("Max. Density = %g 1/m^3\n", maxi);
        printf("Min. Density = %g 1/m^3\n", mini);
    }

    // Here we save at each step if this option has been choosed
    if(g_config->save_step_output) {
        SaveOutputFiles(g_config->output_format, iteration);
        printf("Output number %d has been saved\n", iteration);
    }

    if(fabs(g_config->time - g_config->tf) / fabs(g_config->tf) < SMALL) {
        // Compute the various currents on the various defined contacts
        Compute_Currents( );
        return 1;
    }

    return 0;
}
