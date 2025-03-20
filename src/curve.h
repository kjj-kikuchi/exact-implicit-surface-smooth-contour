//
//  curve.h
//  piecewise-quadraic-implicit-surface
//
//  Created by 菊池祐作 on 2025/02/12.
//

#ifndef curve_h
#define curve_h

#include <Eigen/Core>

struct Curve
{
    std::vector<Eigen::Vector3d> P;
    std::vector<std::vector<Eigen::Vector2i>> L;
};

#endif /* curve_h */
