/*
 * @file mesh_align.h
 * @author Fei Zhu
 * @brief align boundary of mesh B to mesh A according to signed distance
 */

#ifndef MESH_ALIGN_H_
#define MESH_ALIGN_H_

#include "Physika_Core/Vectors/vector_3d.h"
#include "Physika_Core/Arrays/array_Nd.h"
#include "Physika_Geometry/Cartesian_Grids/grid.h"
#include "Physika_Geometry/Boundary_Meshes/surface_mesh.h"

//utility method: compute signed distance to the input mesh
void makeLevelSet(const Physika::SurfaceMesh<float> &surface_mesh, const Physika::Grid<float,3> &grid,
                  Physika::ArrayND<float,3> &signed_dist, unsigned int exact_band = 1);
//generate gradient field from given signed distance field
void makeLevelSetGradient(const Physika::Grid<float,3> &grid, const Physika::ArrayND<float,3> &signed_dist, Physika::ArrayND<Physika::Vector<float,3>,3> &signed_dist_grad);

//adjust the boundary of source_mesh so that it conforms to that of target_mesh
void alignMesh(Physika::SurfaceMesh<float> &target_mesh, Physika::SurfaceMesh<float> &source_mesh, unsigned int level_set_res = 100);

#endif //MESH_ALIGN_H_
