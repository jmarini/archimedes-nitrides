#include "material.h"
#include "mesh.h"


Material material_node(int i, int j) {
    return g_materials[g_mesh->nodes[i][j].material_id];
}


char* mc_material_name(Material *material) {
    switch(material->id) {
        case SILICON: return "Si";
        case GAAS: return "GaAs";
        case GERMANIUM: return "Ge";
        case INSB: return "InSb";
        case ALSB: return "AlSb";
        case ALXINXSB: return "AlxInxSb";
        case ALXIN1XSB: return "AlxIn1xSb";
        case ALAS: return "AlAs";
        case ALP: return "AlP";
        case GAP: return "GaP";
        case GASB: return "GaSb";
        case INAS: return "InAs";
        case INP: return "InP";
        case INXGA1XAS: return "InxGa1xAs";
        case INXAL1XAS: return "InxAl1xAs";
        case INXGAXXAS: return "InxGaXxAs";
        case GAN: return "GaN";
        default: return "Unknown Material";
    }
}
