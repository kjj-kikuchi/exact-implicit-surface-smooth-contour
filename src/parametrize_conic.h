//
//  parametrize_conic.h
//  piecewise-quadraic-implicit-surface
//
//  Created by 菊池祐作 on 2025/03/06.
//

#ifndef parametrize_conic_h
#define parametrize_conic_h

#include <Eigen/Dense>
#include "polynomial/bivariate_quadratic.h"
#include "polynomial/bivariate_linear.h"
#include "curve.h"

double const ZERO = 1.0e-5;
int const NUM_SEGMENTS = 5;
double const INF = std::numeric_limits<double>::infinity();


inline bool matrix_is_singular(Eigen::Vector2d const& sigma)
{
    return fabs( sigma(1) / sigma(0) ) <= ZERO;
}

inline bool equal_zero(double const& value)
{ return fabs(value) <= ZERO; }

inline bool equal_zero(double const& value, Eigen::Vector2d const& sigma, Eigen::Vector2d const& b)
{
    int max_idx;
    Eigen::Vector4d coeffs(sigma(0), sigma(1), b(0), b(1));
    coeffs.cwiseAbs().maxCoeff(&max_idx);
    return fabs(value / coeffs(max_idx)) <= ZERO;
}

inline bool equal_zero(double const& value, Eigen::Vector2d const& sigma, Eigen::Vector2d const& b, double const& c)
{
    int max_idx;
    Eigen::Vector<double, 5> coeffs(sigma(0), sigma(1), b(0), b(1), c);
    coeffs.cwiseAbs().maxCoeff(&max_idx);
    return fabs(value / coeffs(max_idx)) <= ZERO;
}

inline Eigen::Vector3d lifting(Eigen::Vector2d const& v2,
                               BivariateLinear const& removed_term,
                               int const& rt_idx)
{
    Eigen::Vector3d result;
    double lifting_coodinate = removed_term.coeffs(0) * v2(0) + removed_term.coeffs(1) * v2(1) + removed_term.coeffs(2);
    if (rt_idx == 2) {
        result(0) = v2(0);
        result(1) = v2(1);
        result(2) = lifting_coodinate;
    }
    if (rt_idx == 1) {
        result(0) = v2(1);
        result(1) = lifting_coodinate;
        result(2) = v2(0);
    }
    if (rt_idx == 0) {
        result(0) = lifting_coodinate;
        result(1) = v2(0);
        result(2) = v2(1);
    }

    return result;
}


enum class ConicType
{
    ellipse,
    hyperbola_1,
    hyperbola_2,
    intersecting_lines,
    parabola,
    parallel_lines,
    invalid
};


struct ParametrizeConic
{
    ConicType conic_type;
    Eigen::Matrix2d A;
    Eigen::Vector2d b;
    double c;
    Eigen::Matrix2d U;
    Eigen::Vector2d sigma;
    double c_hat;
    Eigen::Vector2d k;
    Eigen::Vector2d b_hat;
    
    ParametrizeConic()
    {
        c_hat = -INF;
        k.setConstant(-INF);
        b_hat.setConstant(-INF);
        conic_type = ConicType::invalid;
    }
    // 2次曲線を識別
    void identify_conic(BivariateQuadratic const& algebraic_curve)
    {
        // Convert into matrix form
        A << 2 * algebraic_curve.coeffs(0),
        algebraic_curve.coeffs(2),
        algebraic_curve.coeffs(2),
        2 * algebraic_curve.coeffs(1);
        b << algebraic_curve.coeffs(3), algebraic_curve.coeffs(4);
        c = algebraic_curve.coeffs(5);

        // Diagonalize the equation
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> es(A);
        U = es.eigenvectors().transpose();
        sigma = es.eigenvalues();

        if (fabs(sigma(0)) < fabs(sigma(1)))    // Let eigen value be |σ1| > |σ2|
        {
            Eigen::Matrix2d P{{0, 1}, {1, 0}};
            sigma = P * sigma;
            U = P * U;
        }

        // Identyfy a conic type
        // A is non-singular
        if (! matrix_is_singular(sigma)) {
            c_hat = c - (b.transpose() * A.inverse() * b)(0, 0) / 2.0;
            // ensure c_hat >= 0 by negating eigen value
            if (c_hat < 0) {
                c_hat = -c_hat;
                sigma = -sigma;
            }

            if (equal_zero(c_hat, sigma, b) && sigma(0) * sigma(1) < 0) conic_type = ConicType::intersecting_lines;
            else if (c_hat > 0) {
                k(0) = std::sqrt( 2 * c_hat / std::fabs(sigma(0)) );
                k(1) = std::sqrt( 2 * c_hat / std::fabs(sigma(1)) );

                if (sigma(0) < 0 && sigma(1) < 0) conic_type = ConicType::ellipse;
                else if (sigma(0) * sigma(1) < 0) {
                    if (sigma(0) * c_hat < 0) conic_type = ConicType::hyperbola_1;
                    else if (sigma(1) * c_hat < 0) conic_type = ConicType::hyperbola_2;
                }
            }
            else conic_type = ConicType::invalid;
        } else {     // A is singular
            b_hat = U * b;
            // negate terms
            if (sigma(0) < 0) {
                sigma = -sigma;
                b_hat = -b_hat;
                c = -c;
            }

            if (!equal_zero(b(1), sigma, b_hat, c)) conic_type = ConicType::parabola;
            else if (b_hat(0)*b_hat(0) > 2*sigma(0)*c) conic_type = ConicType::parallel_lines;
            else conic_type = ConicType::invalid;
        }
    }


    // 有理パラメトリック曲線
    Eigen::Vector2d rational_parametric_curve(double const t, int const j)
    {
        Eigen::Vector2d r;

        if (conic_type == ConicType::ellipse) {
            r(0) = k(0) * (1 - t*t) / (1 + t*t);
            r(1) = k(1) * 2 * t / (1 + t*t);
            if (j == 1) r = -r;
        }
        else if (conic_type == ConicType::hyperbola_1) {
            r(0) = k(0) * (1 + t*t) / (t*t - 1);
            r(1) = k(1) * 2 * t / (t*t - 1);
            if (j == 1) r = -r;
        }
        else if (conic_type == ConicType::hyperbola_2) {
            r(0) = k(0) * 2 * t / (t*t - 1);
            r(1) = k(1) * (1 + t*t) / (t*t - 1);
            if (j == 1) r = -r;
        }
        else if (conic_type == ConicType::intersecting_lines) {
            r(0) = t;
            r(1) = t * std::sqrt(- sigma(0) / sigma(1));
            if (j == 1) r(1) = -r(1);
        }
        else if (conic_type == ConicType::parabola) {
            r(0) = t;
            r(1) = ( - sigma(0) * t*t / 2.0 - b_hat(0) * t - c) / b_hat(1);
        }
        else if (conic_type == ConicType::parallel_lines) {
            if (j == 0) r(0) = ( - b_hat(0) + std::sqrt( b_hat(0)*b_hat(0) - 2*sigma(0)*c ) ) / sigma(0);
            else r(0) = ( - b_hat(0) - std::sqrt( b_hat(0)*b_hat(0) - 2*sigma(0)*c ) ) / sigma(0);
            r(1) = t;
        }

        return r;
    }


    // サンプリング
    void parametrize_conic(BivariateLinear const& removed_term,
                           int const& rt_idx,
                           std::array<std::vector<Eigen::Vector2d>, 2> const& interval,
                           Curve& curve)
    {
        Curve contour_segment;
        for (int j = 0; j < 2; j++) {
            std::vector<Eigen::Vector2d> interval_j = interval[j];

            for (auto const& itvl : interval_j) {
                double spacing = (itvl(1) - itvl(0)) / (double)NUM_SEGMENTS;
                std::vector<Eigen::Vector2i> line_segment;

                for (int i = 0; i <= NUM_SEGMENTS; i++) {
                    double t = itvl(0) + i * spacing;
                    // 2次曲線のパラメータ化
                    Eigen::Vector2d r = rational_parametric_curve(t, j);
                    // 座標変換
                    Eigen::Vector2d x2;
                    if (conic_type == ConicType::ellipse ||
                        conic_type == ConicType::hyperbola_1 || conic_type == ConicType::hyperbola_2 ||
                        conic_type == ConicType::intersecting_lines) {
                        x2 = U.transpose() * r - A.inverse() * b;
                    }
                    else if (conic_type == ConicType::parabola || conic_type == ConicType::parallel_lines) {
                        x2 = U.transpose() * r;
                    }
                    // リフティング
                    Eigen::Vector3d x3 = lifting(x2, removed_term, rt_idx);
                    // 曲線を追加
                    curve.P.push_back(x3);
                    if (i != 0) line_segment.push_back({curve.P.size() - 2, curve.P.size() - 1});
                }
                curve.L.push_back(line_segment);
            }
        }
    }

};


#endif /* parametrize_conic_h */
