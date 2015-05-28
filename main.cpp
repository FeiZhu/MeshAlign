/*
 * @file main.cpp
 * @author Fei Zhu
 * @brief test mesh align alogrithm
 */

#include <string>
#include "Physika_Geometry/Boundary_Meshes/surface_mesh.h"
#include "Physika_IO/Surface_Mesh_IO/surface_mesh_io.h"
#include "mesh_align.h"
using namespace std;

int main()
{
    string target_mesh_file("plant_fine.obj");
    string source_mesh_file("plant_coarse.obj");
    Physika::SurfaceMesh<float> target_mesh, source_mesh;
    Physika::SurfaceMeshIO<float>::load(target_mesh_file,&target_mesh);
    Physika::SurfaceMeshIO<float>::load(source_mesh_file,&source_mesh);
    alignMesh(target_mesh,source_mesh);
    Physika::SurfaceMeshIO<float>::save(target_mesh_file,&target_mesh);
    Physika::SurfaceMeshIO<float>::save(source_mesh_file,&source_mesh);
    return 0;
}
