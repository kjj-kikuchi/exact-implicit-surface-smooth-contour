//
//  bcc_lattice.h
//  piecewise-quadraic-implicit-surface
//
//  Created by 菊池祐作 on 2025/03/06.
//

#ifndef bcc_lattice_h
#define bcc_lattice_h

#include <iostream>
#include <vector>
#include <map>
#include <Eigen/Core>
#include <igl/list_to_matrix.h>

// 格子点の四面体を追加
inline void add_tetrahedras_for_dualpoint(int nx, int ny, int nz,
                                          int ix, int iy, int iz,
                                          std::vector<Eigen::Vector4i>& tet)
{
    int base = ix + (nx * iy) + (nx * ny * iz);

    // primal grid point のインデックス
    int p1 = base + 1;                   // (ix+1,  iy,     iz)
    int p2 = base + nx;                  // (ix,    iy+1,   iz)
    int p3 = base + nx * ny;             // (ix,    iy,     iz+1)
    int p4 = base + nx + 1;              // (ix+1,  iy+1,   iz)
    int p5 = base + nx * ny + nx;        // (ix,    iy+1,   iz+1)
    int p6 = base + nx * ny + 1;         // (ix+1,  iy,     iz+1)
    int p7 = base + nx * ny + nx + 1;    // (ix+1,  iy+1,   iz+1)
    
    // dual grid point のインデックス
    int d0 = base + (nx * ny * nz);
    int d1 = base + 1 + (nx * ny * nz);
    int d2 = base + nx + (nx * ny * nz);
    int d3 = base + nx * ny + (nx * ny * nz);

    // 四面体を追加
    tet.push_back(Eigen::Vector4i(d0, p1, p4, d1));
    tet.push_back(Eigen::Vector4i(d0, p4, p7, d1));
    tet.push_back(Eigen::Vector4i(d0, p7, p6, d1));
    tet.push_back(Eigen::Vector4i(d0, p6, p1, d1));

    tet.push_back(Eigen::Vector4i(d0, p2, p5, d2));
    tet.push_back(Eigen::Vector4i(d0, p5, p7, d2));
    tet.push_back(Eigen::Vector4i(d0, p7, p4, d2));
    tet.push_back(Eigen::Vector4i(d0, p4, p2, d2));

    tet.push_back(Eigen::Vector4i(d0, p3, p6, d3));
    tet.push_back(Eigen::Vector4i(d0, p6, p7, d3));
    tet.push_back(Eigen::Vector4i(d0, p7, p5, d3));
    tet.push_back(Eigen::Vector4i(d0, p5, p3, d3));
}

// BCC lattice を作成
inline void construct_lattice(Eigen::MatrixXd const& V,
                              double const& cell_resolution,
                              Eigen::Vector3i& n,
                              Eigen::MatrixXd& P,
                              Eigen::MatrixXi& Tet)
{
    // bounding box を計算
    Eigen::Vector3d min = V.colwise().minCoeff();
    Eigen::Vector3d max = V.colwise().maxCoeff();

    // 格子サイズを計算
    //    cell_size *= std::max( {max(0) - min(0), max(1) - min(1), max(2) - min(2)} );
    Eigen::Vector3d cell_size = Eigen::Vector3d::Ones() * cell_resolution;
    cell_size = cell_size.cwiseProduct(max - min);
    n(0) = std::ceil((max(0) - min(0)) / cell_size(0)) + 2;
    n(1) = std::ceil((max(1) - min(1)) / cell_size(1)) + 2;
    n(2) = std::ceil((max(2) - min(2)) / cell_size(2)) + 2;

    // grid points を追加
    P.resize(n(0) * n(1) * n(2) * 2, 3);
    int point_idx = 0;
    for (int iz = 0; iz < n(2); iz++) {
        for (int iy = 0; iy < n(1); iy++) {
            for (int ix = 0; ix < n(0); ix++) {
                P.row(point_idx) = min + Eigen::Vector3d(ix-0.5, iy-0.5, iz-0.5).cwiseProduct(cell_size);
                P.row(point_idx + n(0) * n(1) * n(2)) = min + Eigen::Vector3d(ix, iy, iz).cwiseProduct(cell_size);
                point_idx ++;
            }
        }
    }

    // 四面体グリッドを作成
    std::vector<Eigen::Vector4i> tet_list;
    tet_list.reserve( (n(0) * n(1) * n(2)) * (3 * 4));
    for (int iz = 0; iz < n(2) - 1; iz++) {
        for (int iy = 0; iy < n(1) - 1; iy++) {
            for (int ix = 0; ix < n(0) - 1; ix++) {
                add_tetrahedras_for_dualpoint(n(0), n(1), n(2), ix, iy, iz, tet_list);
            }
        }
    }
    Tet.resize(tet_list.size(), 4);
    for (int i = 0; i < tet_list.size(); i++) {
        Tet.row(i) = tet_list[i];
    }
    //    igl::list_to_matrix(tet_list, Tet);
}


inline void enumetate_mid_points(Eigen::MatrixXd const& P,
                                 Eigen::MatrixXi const& Tet,
                                 Eigen::MatrixXd& midP,
                                 std::map<std::pair<int, int>, int>& midP_idx,
                                 std::vector<Eigen::Vector<int, 10>>& control_idx)
{
    std::vector<Eigen::Vector3d> mid_point_list;
    control_idx.resize(Tet.rows());
    for (int i = 0; i < Tet.rows(); i++) {
        Eigen::Vector<int, 10> indices;     // 四面体の10個の制御点のインデックス
        indices.head(4) = Tet.row(i);

        int edges[6][2] = {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};
        for (int e = 0; e < 6; e++) {
            int v0 = Tet(i, edges[e][0]);
            int v1 = Tet(i, edges[e][1]);

            if (v0 > v1) std::swap(v0, v1); // 常に小さいインデックスを先に

            auto edge_key = std::make_pair(v0, v1);
            if ( ! midP_idx.contains(edge_key) ) {
                midP_idx[edge_key] = (int)mid_point_list.size();
                mid_point_list.push_back( (P.row(v0) + P.row(v1)) / 2.0 );
            }

            indices[4 + e] = (int)P.rows() + midP_idx[edge_key];
        }

        control_idx[i] = indices;
    }

    
    midP.resize(mid_point_list.size(), 3);
    for (int i = 0; i < mid_point_list.size(); i++) {
        midP.row(i) = mid_point_list[i];
    }
}


#endif /* bcc_lattice_h */
