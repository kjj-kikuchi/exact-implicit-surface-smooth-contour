//
//  quadratic_implicit_surface.h
//  piecewise-quadraic-implicit-surface
//
//  Created by 菊池祐作 on 2025/03/06.
//

#ifndef quadratic_implicit_surface_h
#define quadratic_implicit_surface_h

#include <Eigen/Core>
#include <igl/signed_distance.h>
#include "polynomial/trivariate_quadratic.h"


// メッシュの区分二次陰関数曲面を計算 (四面体内の陰関数の係数を計算）
inline void compute_implicit_surface(Eigen::MatrixXd const& P,
                                     Eigen::VectorXd const& S,
                                     Eigen::MatrixXi const& Tet,
                                     std::vector<Eigen::Vector<int, 10>> const& control_idx,
                                     std::vector<TrivariateQuadratic>& implicit_surfaces)
{
    // 四面体内の陰関数の係数を計算
    implicit_surfaces.resize(Tet.rows());
    for (int i = 0; i < Tet.rows(); i++)
    {
        Eigen::Matrix<double, 10, 10> Mat;
        Eigen::Vector<double, 10> func_value;
        for (int j = 0; j < 4; j++)
        {
            Eigen::Vector3d x = P.row(Tet(i, j));
            Mat.row(j) <<
            x(0)*x(0), x(1)*x(1), x(2)*x(2),
            x(0)*x(1), x(1)*x(2), x(2)*x(0),
            x(0), x(1), x(2), 1;

            func_value(j) = S(Tet(i, j));
        }

        int edges[6][2] = {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};
        for (int e = 0; e < 6; e++)
        {
            Eigen::Vector3d x = ( P.row(Tet(i, edges[e][0])) + P.row(Tet(i, edges[e][1])) ) / 2.0;
            Mat.row(4 + e) <<
            x(0)*x(0), x(1)*x(1), x(2)*x(2),
            x(0)*x(1), x(1)*x(2), x(2)*x(0),
            x(0), x(1), x(2), 1;

            func_value(4 + e) = S(control_idx[i](4 + e));
        }
        implicit_surfaces[i].coeffs = Mat.inverse() * func_value;
    }
}


// ベジェ係数を計算
void compute_bezier_coefficients(Eigen::MatrixXi const& Tet,
                                 Eigen::VectorXd const& GS,
                                 Eigen::VectorXd const& midS,
                                 std::map<std::pair<int, int>, int> const& midP_idx,
                                 std::vector<Eigen::Vector<double, 10>>& bezier_coeffs)
{
    bezier_coeffs.resize(Tet.rows());
    // BCC lattice の格子点 (0 ~ P.rows()-1)
    for (int i = 0; i < Tet.rows(); i++)
    {
        for (int j = 0; j < 4; j++) bezier_coeffs[i](j) = GS(Tet(i, j));

        int edges[6][2] = {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};
        for (int e = 0; e < 6; e++)
        {
            int v0 = Tet(i, edges[e][0]);
            int v1 = Tet(i, edges[e][1]);
            if (v0 > v1) std::swap(v0, v1);
            auto edge = std::make_pair(v0, v1);
            int edge_idx = midP_idx.at(edge);
            bezier_coeffs[i](4 + e) = 2 * midS(edge_idx) - (GS(edge.first) + GS(edge.second)) / 2.0;
        }
    }
}

void compute_bezier_coefficients(Eigen::MatrixXd const& P,
                                 Eigen::MatrixXd const& midP,
                                 Eigen::MatrixXi const& Tet,
                                 std::vector<Eigen::Vector<int, 10>> const& control_idx,
                                 TrivariateQuadratic const& implicit_surface,
                                 std::map<std::pair<int, int>, int> const& midP_idx,
                                 std::vector<Eigen::Vector<double, 10>>& bezier_coeffs)
{
    Eigen::VectorXd GS(P.rows());
    Eigen::VectorXd midS(midP.rows());

    Eigen::Vector<double, 10> f = implicit_surface.coeffs;
    for (int i = 0; i < GS.size(); i++)
    {
        GS(i) =
        f(0) * P(i, 0) * P(i, 0) +
        f(1) * P(i, 1) * P(i, 1) +
        f(2) * P(i, 2) * P(i, 2) +
        f(3) * P(i, 0) * P(i, 1) +
        f(4) * P(i, 1) * P(i, 2) +
        f(5) * P(i, 2) * P(i, 0) +
        f(6) * P(i, 0) +
        f(7) * P(i, 1) +
        f(8) * P(i, 2) +
        f(9);
    }
    for (int i = 0; i < midS.size(); i++)
    {
        midS(i) =
        f(0) * midP(i, 0) * midP(i, 0) +
        f(1) * midP(i, 1) * midP(i, 1) +
        f(2) * midP(i, 2) * midP(i, 2) +
        f(3) * midP(i, 0) * midP(i, 1) +
        f(4) * midP(i, 1) * midP(i, 2) +
        f(5) * midP(i, 2) * midP(i, 0) +
        f(6) * midP(i, 0) +
        f(7) * midP(i, 1) +
        f(8) * midP(i, 2) +
        f(9);
    }
    bezier_coeffs.resize(Tet.rows());
    // BCC lattice の格子点 (0 ~ P.rows()-1)
    for (int i = 0; i < Tet.rows(); i++)
    {
        Eigen::Vector<double, 10> f;
        for (int j = 0; j < 4; j++) bezier_coeffs[i](j) = GS(Tet(i, j));

        int edges[6][2] = {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};
        for (int e = 0; e < 6; e++)
        {
            int v0 = Tet(i, edges[e][0]);
            int v1 = Tet(i, edges[e][1]);
            if (v0 > v1) std::swap(v0, v1);
            auto edge = std::make_pair(v0, v1);
            int edge_idx = midP_idx.at(edge);
            bezier_coeffs[i](4 + e) = 2 * midS(edge_idx) - (GS(edge.first) + GS(edge.second)) / 2.0;
        }
    }
}

void compute_bezier_coefficients(Eigen::VectorXd const& GS,
                                 Eigen::VectorXd const& midS,
                                 std::map<std::pair<int, int>, int> const& midP_idx,
                                 Eigen::VectorXd& bezier_coeffs)
{
    bezier_coeffs.resize(GS.rows() + midS.rows());
    // BCC lattice の格子点 (0 ~ P.rows()-1)
    for (int i = 0; i < GS.rows(); i++) bezier_coeffs(i) = GS(i);
    // 格子点の中点 (P.rows() ~ P.rows() + midP.rows()-1)
    for (auto [edge, idx] : midP_idx)
    {
        bezier_coeffs(GS.size() + idx) = 2 * midS(idx) - (GS(edge.first) + GS(edge.second)) / 2.0;
    }
}

// signed distanceを計算
inline void get_signed_distance(Eigen::MatrixXd const& P,
                                Eigen::MatrixXd const& V,
                                Eigen::MatrixXi const& F,
                                Eigen::VectorXd& S)
{
    Eigen::VectorXi I;
    Eigen::MatrixXd C;
    Eigen::MatrixXd N;
    igl::signed_distance(P, V, F, igl::SIGNED_DISTANCE_TYPE_FAST_WINDING_NUMBER, S, I, C, N);
}
#endif /* quadratic_implicit_surface_h */
