//
//  trivariate_quadratic.h
//  piecewise-quadraic-implicit-surface
//
//  Created by 菊池祐作 on 2025/03/06.
//

#ifndef trivariate_quadratic_h
#define trivariate_quadratic_h

#include "trivariate_linear.h"

struct TrivariateQuadratic
{
    Eigen::Vector<double, 10> coeffs;

    TrivariateLinear derivative(int axis) const
    {
        TrivariateLinear result;
        result.coeffs(axis) = 2 * coeffs(axis);
        result.coeffs((axis + 1)%3) = coeffs(axis + 3);
        result.coeffs((axis + 2)%3) = coeffs((axis + 2)%3 + 3);
        result.coeffs(3) = coeffs(axis + 6);

        return result;
    }
};

#endif /* trivariate_quadratic_h */
