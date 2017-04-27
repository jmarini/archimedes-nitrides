#include "mesh.h"

#include <stdio.h>

#include "random.h"


int mc_build_mesh(Mesh *mesh) {
    // saving into local variables for readability
    int nx = mesh->nx,
        ny = mesh->ny;
    double dx = mesh->dx,
           dy = mesh->dy;

    // definition of mesh node coordinates
    mesh->num_nodes = (nx + 1) * (ny + 1);
    for(int i = 0; i < nx + 1; ++i) {
        for(int j = 0; j < ny + 1; ++j) {
            int index = j * (nx + 1) + i;

            mesh->coordinates[index] = (Vec2){.x=(double)(i) * dx,
                                              .y=(double)(j) * dy};
            mesh->nodes[i+1][j+1].index = (Index){.i=i+1, .j=j+1};
        }
    }

    // definition of mesh triangles
    mesh->num_triangles = 2 * nx * ny;
    for(int i = 0; i < nx; ++i) {
        for(int j = 0; j < ny; ++j) {
            // bottom up triangle
            int index = j * nx + i;
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


int mc_save_mesh(Mesh *mesh, char *filename) {
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
        fprintf(fp, "%g %g 0\n", mesh->coordinates[i].x, mesh->coordinates[i].y);
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


int mc_is_boundary_insulator(int direction, int index) {
    return g_mesh->edges[direction][index].boundary == boundary_t.INSULATOR;
}


int mc_is_boundary_schottky(int direction, int index) {
    return g_mesh->edges[direction][index].boundary == boundary_t.SCHOTTKY;
}


int mc_is_boundary_ohmic(int direction, int index) {
    return g_mesh->edges[direction][index].boundary == boundary_t.OHMIC;
}


int mc_is_boundary_vacuum(int direction, int index) {
    return g_mesh->edges[direction][index].boundary == boundary_t.VACUUM;
}


int mc_is_boundary_contact(int direction, int index) {
    return mc_is_boundary_schottky(direction, index)
        || mc_is_boundary_ohmic(direction, index);
}


Node * mc_node(int i, int j)  { return &(g_mesh->nodes[i][j]); }


Node * mc_node_s(Index index) { return &(g_mesh->nodes[index.i][index.j]); }


Node * mc_get_particle_node(Particle *p) {
    return mc_node_s(mc_particle_coords(p));
}


Vec2 mc_random_location_in_node(Node *node) {
    double x = g_mesh->dx * ((double)node->i - rnd());
    double y = g_mesh->dy * ((double)node->j - rnd());

    return (Vec2){.x=x, .y=y};
}
