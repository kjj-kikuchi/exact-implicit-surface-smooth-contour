//
//  trimming_conic.h
//  piecewise-quadraic-implicit-surface
//
//  Created by 菊池祐作 on 2025/03/06.
//

#ifndef trimming_conic_h
#define trimming_conic_h

#include <iostream>
#include <Eigen/LU>
#include "parametrize_conic.h"
#include "polynomial/trivariate_linear.h"
#include "polynomial/bivariate_quadratic.h"
#include "polynomial/bivariate_linear.h"
#include "polynomial/rational_parametric_curve.h"
#include "curve.h"

// 2つの平面の範囲の積集合を求める
inline std::vector<Eigen::Vector2d> intersect_two_sets(std::vector<Eigen::Vector2d> const& set1,
                                                       std::vector<Eigen::Vector2d> const& set2)
{
    std::vector<Eigen::Vector2d> result;

    for (auto const& itvl1 : set1) {
        for (auto const& itvl2 : set2) {
            double lower = std::max(itvl1(0), itvl2(0));
            double upper = std::min(itvl1(1), itvl2(1));
            if (lower < upper) result.push_back({lower, upper});
        }
    }
    return result;
}

// 複数の平面の範囲の積集合を求める
inline std::vector<Eigen::Vector2d> get_intersection_set(std::vector<std::vector<Eigen::Vector2d>> const& plane_intervals)
{
    if (plane_intervals.empty()) return {};

    // 初期範囲（全空間）
    std::vector<Eigen::Vector2d> intersection = { {-INF, INF} };

    // 各平面の範囲と順番に積集合を計算
    for (auto const& plane : plane_intervals)
    {
        intersection = intersect_two_sets(intersection, plane);
        if (intersection.empty()) return {}; // 積集合が空になったら終了
    }

    return intersection;
}

inline bool coeff_equal_zero(int const& idx, Eigen::VectorXd const& coeffs)
{
    int max_idx;
    coeffs.cwiseAbs().maxCoeff(&max_idx);
    return fabs(coeffs(idx) / coeffs.mean()) <= ZERO;
}


// トリミング
inline std::array<std::vector<Eigen::Vector2d>, 2>
trimming_conic(std::vector<Eigen::Vector3d> const& point,
               std::vector<Eigen::Vector3d> const& normal,
               ParametrizeConic const& pc,
               BivariateLinear const& removed_term,
               int const& rt_idx)
{
    std::array<std::vector<Eigen::Vector2d>, 2> result;
    // 曲線を二つに分けて計算  j == 1 のとき符号反転
    for (int j = 0; j < 2; j++) {
        if (j == 1 && pc.conic_type == ConicType::parabola) break;
        std::vector<std::vector<Eigen::Vector2d>> plane_intervals;     // plane_intervals[i] : 平面 i の範囲

        // 四面体の各平面に対して範囲を計算
        for (int i = 0; i < point.size(); i++) {
            Eigen::Vector3d p = point[i];
            Eigen::Vector3d n = normal[i];

            // 有理パラメトリック曲線
            RationalParametricCurve2 r;
            if (pc.conic_type == ConicType::ellipse) r.ellipse(pc.k);
            if (pc.conic_type == ConicType::hyperbola_1) r.hyperbola_1(pc.k);
            if (pc.conic_type == ConicType::hyperbola_2) r.hyperbola_2(pc.k);
            if (pc.conic_type == ConicType::intersecting_lines) r.intersecting_lines(pc.sigma);
            if (pc.conic_type == ConicType::parabola) r.parabola(pc.sigma(0), pc.b_hat, pc.c);
            if (pc.conic_type == ConicType::parallel_lines) r.parallel_lines(pc.sigma(0), pc.b_hat(0), pc.c, j);

            if (j == 1 && (pc.conic_type == ConicType::ellipse ||
                           pc.conic_type == ConicType::hyperbola_1 || pc.conic_type == ConicType::hyperbola_2)) {
                r.numer_x.coeffs = - r.numer_x.coeffs;
                r.numer_y.coeffs = - r.numer_y.coeffs;
            }
            else if (j == 1 && pc.conic_type == ConicType::intersecting_lines) {
                r.numer_y.coeffs = - r.numer_y.coeffs;
            }

            // デカルト座標空間に変換
            RationalParametricCurve2 x2;
            if (pc.conic_type == ConicType::ellipse ||
                pc.conic_type == ConicType::hyperbola_1 || pc.conic_type == ConicType::hyperbola_2 ||
                pc.conic_type == ConicType::intersecting_lines) {
                x2 = pc.U.transpose() * r - (pc.A.inverse() * pc.b) * r.denom;
            }
            else if (pc.conic_type == ConicType::parabola || pc.conic_type == ConicType::parallel_lines) {
                x2 = pc.U.transpose() * r;
            }

            // 2次元から3次元にリフティング
            RationalParametricCurve3 x3 = lifting(x2, removed_term, rt_idx);

            // 範囲tの2次方程式    n * (x_d(t) - x_n(t) p) <= 0
            UnivariateQuadratic dot;
            dot.coeffs =
            n(0) * (x3.numer_x.coeffs - x3.denom.coeffs * p(0)) +
            n(1) * (x3.numer_y.coeffs - x3.denom.coeffs * p(1)) +
            n(2) * (x3.numer_z.coeffs - x3.denom.coeffs * p(2));

            if (pc.conic_type == ConicType::hyperbola_1 ||
                pc.conic_type == ConicType::hyperbola_2) dot.coeffs = - dot.coeffs;   // 双曲線の場合、不等号を反転

            // 区間を計算
            std::vector<Eigen::Vector2d> interval;
            if (!equal_zero(dot.coeffs(0))) { // 2次の場合
                double D = dot.coeffs(1) * dot.coeffs(1) - 4 * dot.coeffs(0) * dot.coeffs(2);   // 判定式

                if (D >= 0) {   // 解が存在
                    double t1 = (-dot.coeffs(1) + std::sqrt(D)) / (2 * dot.coeffs(0));
                    double t2 = (-dot.coeffs(1) - std::sqrt(D)) / (2 * dot.coeffs(0));
                    //double t2 = dot.coeffs(2) / (dot.coeffs(0) * t1);

                    if (dot.coeffs(0) > 0) interval.push_back({t2, t1});
                    else if (dot.coeffs(0) < 0) {
                        interval.push_back({-INF, t1});
                        interval.push_back({t2, INF});
                    }
                    plane_intervals.push_back(interval);
                } else { // 解が存在しない
                    if (dot.coeffs(0) > 0) {
                        plane_intervals = {};
                        break;
                    }
                    else if (dot.coeffs(0) < 0) {
                        interval.push_back({-INF, INF});
                        plane_intervals.push_back(interval);
                    }
                }
            } else if (! equal_zero(dot.coeffs(1))) { // 1次の場合
                if (dot.coeffs(1) > 0) interval.push_back({-INF, - dot.coeffs(2) / dot.coeffs(1)});
                else interval.push_back({- dot.coeffs(2) / dot.coeffs(1), INF});
            } else {
                plane_intervals = {};
                break;
            }
            plane_intervals.push_back(interval);
        }

        // 楕円・双曲線のときは t : [-1, 1]とする
        std::vector<Eigen::Vector2d> interval;
        if (pc.conic_type == ConicType::ellipse && plane_intervals.size() != 0) {
            interval.push_back({-1, 1});
            plane_intervals.push_back(interval);
        }
        else if ((pc.conic_type == ConicType::hyperbola_1 || pc.conic_type == ConicType::hyperbola_2) &&
                 plane_intervals.size() != 0) {
            interval.push_back({-0.99999, 0.99999});
            plane_intervals.push_back(interval);
        }

        result[j] = get_intersection_set(plane_intervals);
    }
    return result;
}


#endif /* trimming_conic_h */
