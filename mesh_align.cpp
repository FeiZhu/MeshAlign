/*
 * @file mesh_align.cpp
 * @author Fei Zhu
 * @brief align boundary of mesh B to mesh A according to signed distance
 */

#include <vector>
#include <iostream>
#include "Physika_Core/Utilities/math_utilities.h"
#include "Physika_Core/Range/range.h"
#include "./SDFGen/makelevelset3.h"
#include "mesh_align.h"

void makeLevelSet(const Physika::SurfaceMesh<float> &surface_mesh, const Physika::Grid<float,3> &grid,
                  Physika::ArrayND<float,3> &signed_dist, unsigned int exact_band)
{
    if(surface_mesh.isTriangularMesh() == false)
    {
        std::cout<<"Error: only triangle surface mesh supported!\n";
        return;
    }
    //data converter, generate data needed by SDFGen method
    std::vector<Vec3ui> tri;
    for(unsigned int face_idx = 0; face_idx < surface_mesh.numFaces(); ++face_idx)
    {
        const Physika::SurfaceMeshInternal::Face<float> &face = surface_mesh.face(face_idx);
        Vec3ui triangle;
        for(unsigned int vert_idx = 0; vert_idx < face.numVertices(); ++vert_idx)
        {
            const Physika::BoundaryMeshInternal::Vertex<float> &vertex = face.vertex(vert_idx);
            triangle[vert_idx] = vertex.positionIndex();
        }
        tri.push_back(triangle);
    }
    std::vector<Vec3f> x;
    for(unsigned int vert_idx = 0; vert_idx < surface_mesh.numVertices(); ++vert_idx)
    {
        Physika::Vector<float,3> vert_pos = surface_mesh.vertexPosition(vert_idx);
        Vec3f x_pos;
        for(unsigned int dim = 0; dim < 3; ++dim)
            x_pos[dim] = vert_pos[dim];
        x.push_back(x_pos);
    }
    Vec3f origin;
    Physika::Vector<float,3> grid_min_corner = grid.minCorner();
    for(unsigned int dim  = 0; dim < 3; ++dim)
        origin[dim] = grid_min_corner[dim];
    float dx;
    Physika::Vector<float,3> grid_dx = grid.dX();
    if(!(grid_dx[0] == grid_dx[1] && grid_dx[1] == grid_dx[2]))
    {
        std::cout<<"Error: only support grid with uniform grid size in each direction!\n";
        return;
    }
    dx = grid_dx[0];
    int nx,ny,nz;
    Physika::Vector<unsigned int,3> grid_node_num = grid.nodeNum();
    nx = grid_node_num[0]; ny = grid_node_num[1]; nz = grid_node_num[2];
    Array3f phi;
    make_level_set3(tri,x,origin,dx,nx,ny,nz,phi,exact_band);
    signed_dist.resize(grid_node_num);
    for(unsigned int x = 0; x < nx; ++x)
        for(unsigned int y = 0; y < ny; ++y)
            for(unsigned z = 0; z < nz; ++z)
            {
                Physika::Vector<unsigned int,3> idx(x,y,z);
                signed_dist(idx) = phi(x,y,z);
            }
}

void makeLevelSetGradient(const Physika::Grid<float,3> &grid, const Physika::ArrayND<float,3> &signed_dist, Physika::ArrayND<Physika::Vector<float,3>,3> &signed_dist_grad)
{
    Physika::Vector<unsigned int,3> grid_res = grid.nodeNum();
    Physika::Vector<float,3> grid_dx = grid.dX();
    signed_dist_grad.resize(grid_res);
    for(Physika::Grid<float,3>::NodeIterator iter = grid.nodeBegin(); iter != grid.nodeEnd(); ++iter)
    {
        Physika::Vector<unsigned int,3> node_idx = iter.nodeIndex();
        //compute gradient with central difference
        //for grid node at the boundary, use one side finite difference
        for(unsigned int dim = 0; dim < 3; ++dim)
        {
            unsigned length_scale = 0;
            Physika::Vector<unsigned int,3> pos_idx = node_idx;
            Physika::Vector<unsigned int,3> neg_idx = node_idx;
            if(pos_idx[dim] < grid_res[dim] - 1)
            {
                pos_idx[dim] ++;
                ++length_scale;
            }
            if(neg_idx[dim] > 0)
            {
                neg_idx[dim] --;
                ++length_scale;
            }
            signed_dist_grad(node_idx)[dim] = (signed_dist(pos_idx) - signed_dist(neg_idx))/(length_scale*grid_dx[dim]);
        }
    }
}

void alignMesh(Physika::SurfaceMesh<float> &target_mesh, Physika::SurfaceMesh<float> &source_mesh, unsigned int level_set_res)
{
    Physika::Range<float,3> target_mesh_aabb = target_mesh.axisAlignedBoundingBox();
    Physika::Range<float,3> source_mesh_aabb = source_mesh.axisAlignedBoundingBox();
    //the grid range is the union of the two bounding boxes
    Physika::Range<float,3> grid_range;
    Physika::Vector<float,3> grid_min_corner,grid_max_corner;
    for(unsigned int dim = 0; dim < 3; ++dim)
    {
        grid_min_corner[dim] = Physika::min(target_mesh_aabb.minCorner()[dim],source_mesh_aabb.minCorner()[dim]);
        grid_max_corner[dim] = Physika::max(target_mesh_aabb.maxCorner()[dim],source_mesh_aabb.maxCorner()[dim]);
    }
    grid_range.setMinCorner(grid_min_corner);
    grid_range.setMaxCorner(grid_max_corner);
    Physika::Grid<float,3> grid(grid_range,level_set_res-1);
    Physika::ArrayND<float,3> signed_dist;
    makeLevelSet(target_mesh,grid,signed_dist);
    Physika::ArrayND<Physika::Vector<float,3>,3> signed_dist_grad;
    makeLevelSetGradient(grid,signed_dist,signed_dist_grad);
    //for each vertex of the source mesh, move it in the negative gradient direction
    //the distance to be moved is the distance from the vertex to target mesh
    for(unsigned int vert_idx = 0; vert_idx < source_mesh.numVertices(); ++vert_idx)
    {
        Physika::Vector<float,3> vert_pos = source_mesh.vertexPosition(vert_idx);
        Physika::Vector<unsigned int,3> cell_idx;
        Physika::Vector<float,3> bias_in_cell;
        grid.cellIndexAndBiasInCell(vert_pos,cell_idx,bias_in_cell);
        //the signed distance at this vertex is interpolated from grid
        //gradient at this vertex is also interpolated from grid
        float vert_signed_dist = 0.0f;
        Physika::Vector<float,3> vert_gradient(0.0f);
        for(unsigned int idx_x = 0; idx_x < 2; ++idx_x)
        {
            float weight_x = idx_x == 0 ? 1 - bias_in_cell[0] : bias_in_cell[0];
            for(unsigned int idx_y = 0; idx_y < 2; ++idx_y)
            {
                float weight_y = idx_y == 0 ? 1 - bias_in_cell[1] : bias_in_cell[1];
                for(unsigned int idx_z = 0; idx_z < 2; ++idx_z)
                {
                    Physika::Vector<unsigned int,3> node_idx = cell_idx;
                    node_idx[0] += idx_x; node_idx[1] += idx_y; node_idx[2] += idx_z;
                    float weight_z = idx_z == 0 ? 1 - bias_in_cell[2] : bias_in_cell[2];
                    float weight = weight_x*weight_y*weight_z;
                    vert_signed_dist += weight*signed_dist(node_idx);
                    vert_gradient += weight*signed_dist_grad(node_idx);
                }
            }
        }
        vert_pos += -Physika::abs(vert_signed_dist)*(vert_gradient.normalize());
        //update the vertex position of the source mesh
        source_mesh.setVertexPosition(vert_idx,vert_pos);
    }
}
