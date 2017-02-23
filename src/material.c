#include "material.h"
#include "mesh.h"


Material material_node(int i, int j) {
    return g_materials[g_mesh->nodes[i][j].material];
}
