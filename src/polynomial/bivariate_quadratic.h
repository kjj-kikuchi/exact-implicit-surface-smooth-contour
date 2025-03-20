//
//  bivariate_quadratic.h
//  piecewise-quadraic-implicit-surface
//
//  Created by 菊池祐作 on 2025/03/06.
//

#ifndef bivariate_quadratic_h
#define bivariate_quadratic_h

#include <Eigen/Core>

struct BivariateQuadratic
{
    Eigen::Vector<double, 6> coeffs;
};

#endif /* bivariate_quadratic_h */
