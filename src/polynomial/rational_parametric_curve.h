//
//  rational_parametric_curve.h
//  piecewise-quadraic-implicit-surface
//
//  Created by 菊池祐作 on 2025/03/06.
//

#ifndef rational_parametric_curve_h
#define rational_parametric_curve_h

#include "univariate_quadratic.h"
#include "bivariate_linear.h"

// 3次元空間内の有理パラメトリック曲線
struct RationalParametricCurve3
{
    UnivariateQuadratic numer_x;
    UnivariateQuadratic numer_y;
    UnivariateQuadratic numer_z;
    UnivariateQuadratic denom;
};

// 2次元平面内の有理パラメトリック曲線
struct RationalParametricCurve2
{
    UnivariateQuadratic numer_x;
    UnivariateQuadratic numer_y;
    UnivariateQuadratic denom;

    void ellipse(Eigen::Vector2d const& k)
    {
        numer_x.coeffs = k(0) * Eigen::Vector3d(-1, 0, 1);      // k1 (- t^2 + 1)
        numer_y.coeffs = 2 * k(1) * Eigen::Vector3d(0, 1, 0);   // 2 k2 t
        denom.coeffs = {1, 0, 1};                               // t^2 + 1
    }

    void hyperbola_1(Eigen::Vector2d const& k)
    {
        numer_x.coeffs = k(0) * Eigen::Vector3d(1, 0, 1);       // k1 (t^2 + 1)
        numer_y.coeffs = 2 * k(1) * Eigen::Vector3d(0, 1, 0);   // 2 k2 t
        denom.coeffs = {1, 0, -1};                              // t^2 - 1
    }

    void hyperbola_2(Eigen::Vector2d const& k)
    {
        numer_x.coeffs = 2 * k(0) * Eigen::Vector3d(0, 1, 0);   // 2 k2 t
        numer_y.coeffs = k(1) * Eigen::Vector3d(1, 0, 1);       // k1 (t^2 + 1)
        denom.coeffs = {1, 0, -1};                              // t^2 - 1
    }

    void intersecting_lines(Eigen::Vector2d const& sigma)
    {
        numer_x.coeffs = {0, 1, 0};     // t
        numer_y.coeffs = std::sqrt(- sigma(0) / sigma(1)) * Eigen::Vector3d(0, 1, 0);   // t √(-σ_1/σ_2)
        denom.coeffs = {0, 0, 1};       // 1
    }

    void parabola(double const& sigma1, Eigen::Vector2d const& b_hat, double const& c)
    {
        numer_x.coeffs = {0, 1, 0};     // t
        numer_y.coeffs = Eigen::Vector3d(sigma1, b_hat(0), c) / -b_hat(1);
        denom.coeffs = {0, 0, 1};       // 1
    }

    void parallel_lines(double const& sigma1, double const& b_hat1, double const& c, int const& j)
    {
        if (j == 0)
            numer_x.coeffs = Eigen::Vector3d(0, 0, 1) * (- b_hat1 + std::sqrt(b_hat1*b_hat1 - 2*sigma1*c)) / sigma1;
        else
            numer_x.coeffs = Eigen::Vector3d(0, 0, 1) * (- b_hat1 - std::sqrt(b_hat1*b_hat1 - 2*sigma1*c)) / sigma1;
        numer_y.coeffs = {0, 1, 0};     // t
        denom.coeffs = {0, 0, 1};       // 1
    }
};

// 2Dから3Dに持ち上げる
inline RationalParametricCurve3 lifting(RationalParametricCurve2 const& r2,
                                        BivariateLinear const& removed_term,
                                        int const& rt_idx)
{
    RationalParametricCurve3 result;
    UnivariateQuadratic lifting_coordinate;
    lifting_coordinate.coeffs =
    removed_term.coeffs(0) * r2.numer_x.coeffs +
    removed_term.coeffs(1) * r2.numer_y.coeffs +
    removed_term.coeffs(2) * r2.denom.coeffs;
    if (rt_idx == 2)
    {
        result.numer_x.coeffs = r2.numer_x.coeffs;
        result.numer_y.coeffs = r2.numer_y.coeffs;
        result.numer_z.coeffs = lifting_coordinate.coeffs;
    }
    if (rt_idx == 1)
    {
        result.numer_x.coeffs = r2.numer_y.coeffs;
        result.numer_y.coeffs = lifting_coordinate.coeffs;
        result.numer_z.coeffs = r2.numer_x.coeffs;
    }
    if (rt_idx == 0)
    {
        result.numer_x.coeffs = lifting_coordinate.coeffs;
        result.numer_y.coeffs = r2.numer_x.coeffs;
        result.numer_z.coeffs = r2.numer_y.coeffs;
    }
    result.denom.coeffs = r2.denom.coeffs;

    return result;
}

inline RationalParametricCurve2 operator*(Eigen::Matrix2d const& M,
                                          RationalParametricCurve2 const& r2)
{
    RationalParametricCurve2 result;
    result.numer_x.coeffs = M(0, 0) * r2.numer_x.coeffs + M(0, 1) * r2.numer_y.coeffs;
    result.numer_y.coeffs = M(1, 0) * r2.numer_x.coeffs + M(1, 1) * r2.numer_y.coeffs;
    result.denom.coeffs = r2.denom.coeffs;
    return result;
}

inline RationalParametricCurve2 operator*(Eigen::Vector2d const& vec,
                                          UnivariateQuadratic const& denom)
{
    RationalParametricCurve2 result;
    result.numer_x.coeffs = vec(0) * denom.coeffs;
    result.numer_y.coeffs = vec(1) * denom.coeffs;
    result.denom.coeffs = denom.coeffs;
    return result;
}


inline RationalParametricCurve2 operator-(RationalParametricCurve2 const& r1,
                                          RationalParametricCurve2 const& r2)
{
    RationalParametricCurve2 result;
    result.numer_x.coeffs = r1.numer_x.coeffs - r2.numer_x.coeffs;
    result.numer_y.coeffs = r1.numer_y.coeffs - r2.numer_y.coeffs;
    result.denom.coeffs = r1.denom.coeffs;
    return result;
}

#endif /* rational_parametric_curve_h */
