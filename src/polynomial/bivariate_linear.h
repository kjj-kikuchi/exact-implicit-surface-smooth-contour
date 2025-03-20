//
//  bivariate_linear.h
//  piecewise-quadraic-implicit-surface
//
//  Created by 菊池祐作 on 2025/03/06.
//

#ifndef bivariate_linear_h
#define bivariate_linear_h

#include "bivariate_quadratic.h"

struct BivariateLinear
{
    Eigen::Vector3d coeffs;

    BivariateQuadratic square() const
    {
        BivariateQuadratic result;
        result.coeffs(0) = coeffs(0) * coeffs(0);
        result.coeffs(1) = coeffs(1) * coeffs(1);
        result.coeffs(2) = 2 * coeffs(0) * coeffs(1);
        result.coeffs(3) = 2 * coeffs(0) * coeffs(2);
        result.coeffs(4) = 2 * coeffs(1) * coeffs(2);
        result.coeffs(5) = coeffs(2) * coeffs(2);

        return result;
    }

    BivariateQuadratic multiple_x() const
    {
        BivariateQuadratic result;
        result.coeffs(0) = coeffs(0);
        result.coeffs(1) = 0;
        result.coeffs(2) = coeffs(1);
        result.coeffs(3) = coeffs(2);
        result.coeffs(4) = 0;
        result.coeffs(5) = 0;

        return result;
    }

    BivariateQuadratic multiple_y() const
    {
        BivariateQuadratic result;
        result.coeffs(0) = 0;
        result.coeffs(1) = coeffs(1);
        result.coeffs(2) = coeffs(0);
        result.coeffs(3) = 0;
        result.coeffs(4) = coeffs(2);
        result.coeffs(5) = 0;

        return result;
    }

    BivariateQuadratic convert_to_quadratic() const
    {
        BivariateQuadratic result;
        result.coeffs(0) = 0;
        result.coeffs(1) = 0;
        result.coeffs(2) = 0;
        result.coeffs(3) = coeffs(0);
        result.coeffs(4) = coeffs(1);
        result.coeffs(5) = coeffs(2);

        return result;
    }
};

#endif /* bivariate_linear_h */
