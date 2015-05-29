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
    string input_target_mesh_file("fine.obj");
    string input_source_mesh_file("coarse.obj");
    Physika::SurfaceMesh<float> target_mesh, source_mesh;
    Physika::SurfaceMeshIO<float>::load(input_target_mesh_file,&target_mesh);
    Physika::SurfaceMeshIO<float>::load(input_source_mesh_file,&source_mesh);
    unsigned int grid_res = 100;
    alignMesh(target_mesh,source_mesh,grid_res);
    string output_source_mesh_file("coarse_out.obj");
    Physika::SurfaceMeshIO<float>::save(output_source_mesh_file,&source_mesh);
    return 0;
}
