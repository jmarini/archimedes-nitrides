#include "global_defines.h"
#include "mesh.h"

#include <stdio.h>


int mc_build_mesh(mc_mesh_t *mesh) {
    int i = 0,
        j = 0,
        index = 0;

    // saving into local variables for readability
    int nx = mesh->nx,
        ny = mesh->ny;
    real dx = mesh->dx,
         dy = mesh->dy;

    // definition of mesh node coordinates
    mesh->num_nodes = (nx + 1) * (ny + 1);
    for(i = 0; i < nx + 1; ++i) {
        for(j = 0; j < ny + 1; ++j) {
            index = j * (nx + 1) + i;
            mesh->coordinates[index][0] = (real)(i) * dx;
            mesh->coordinates[index][1] = (real)(j) * dy;

            mesh->info[i+1][j+1].i = i+1;
            mesh->info[i+1][j+1].j = j+1;
        }
    }

    // definition of mesh triangles
    mesh->num_triangles = 2 * nx * ny;
    for(i = 0; i < nx; ++i) {
        for(j = 0; j < ny; ++j) {
            // bottom up triangle
            index = j * nx + i;
            mesh->triangles[index][0] =  j      * (nx + 1) + i - 1;
            mesh->triangles[index][1] =  j      * (nx + 1) + i;
            mesh->triangles[index][2] = (j + 1) * (nx + 1) + i;

            // top down triangle
            index = j * nx + i + nx * ny;
            mesh->triangles[index][0] =  j      * (nx + 1) + i - 1;
            mesh->triangles[index][1] = (j + 1) * (nx + 1) + i;
            mesh->triangles[index][2] = (j + 1) * (nx + 1) + i - 1;
        }
    }

    return 0;
}


int mc_save_mesh(mc_mesh_t *mesh, char *filename) {
    int i = 0;
    FILE *fp;

    fp = fopen(filename, "w");
    if(fp == NULL) {
       printf("Error could not open file '%s'.\n", filename);
       return 1;
    }

    fprintf(fp, "MeshVersionFormatted 1\n\n");
    fprintf(fp, "Dimension\n2\n\n");
    fprintf(fp, "Vertices\n");
    fprintf(fp, "%d\n", mesh->num_nodes);

    for(i = 0; i< mesh->num_nodes; ++i) {
        fprintf(fp, "%g %g 0\n", mesh->coordinates[i][0], mesh->coordinates[i][1]);
    }

    fprintf(fp, "\n");
    fprintf(fp, "Triangles\n");
    fprintf(fp, "%d\n", mesh->num_triangles);
    for(i = 0; i < mesh->num_triangles; ++i) {
        fprintf(fp,"%d %d %d 0\n",
                mesh->triangles[i][0], mesh->triangles[i][1], mesh->triangles[i][2]);
    }

    fclose(fp);

    return 0;
}


mc_node_t * mc_node(int i, int j) {
    return &(g_mesh->info[i][j]);
}


mc_node_t * mc_node_s(Index index) {
    return &(g_mesh->info[index.i][index.j]);
}
