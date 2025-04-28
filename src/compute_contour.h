//
//  compute_contour.h
//  piecewise-quadraic-implicit-surface
//
//  Created by 菊池祐作 on 2025/03/06.
//

#ifndef compute_contour_h
#define compute_contour_h

#include "quadratic_implicit_surface.h"
#include "parametrize_conic.h"
#include "trimming_conic.h"

inline Eigen::Vector3d get_normal(Eigen::Vector3d const& p0, Eigen::Vector3d const& p1, Eigen::Vector3d const& p2)
{ return ((p1-p0).cross(p2-p0)).normalized(); }

inline bool is_trivial_reject_tet(Eigen::Vector<double, 10> const& bezier_coeffs)
{ return (bezier_coeffs.array() > 0).all() || (bezier_coeffs.array() < 0).all(); }

inline void tetrahedral_faces(Eigen::MatrixXd const& P,
                              Eigen::MatrixXi const& Tet,
                              int const& i,
                              std::vector<Eigen::Vector3d>& point,
                              std::vector<Eigen::Vector3d>& normal)
{
    for (int j = 0; j < 4; j++) point[j] = P.row(Tet(i, j));
    normal[0] = get_normal(P.row(Tet(i, 0)), P.row(Tet(i, 1)), P.row(Tet(i, 3)));
    normal[1] = get_normal(P.row(Tet(i, 0)), P.row(Tet(i, 2)), P.row(Tet(i, 1)));
    normal[2] = get_normal(P.row(Tet(i, 0)), P.row(Tet(i, 3)), P.row(Tet(i, 2)));
    normal[3] = get_normal(P.row(Tet(i, 1)), P.row(Tet(i, 2)), P.row(Tet(i, 3)));
}


inline TrivariateLinear orthographic_orthogonality_condition(TrivariateQuadratic const& f,
                                                              Eigen::Vector3d const& camera)
{
    // 勾配ベクトル ∇f(x,y,z)
    TrivariateLinear partial_x = f.derivative(0);
    TrivariateLinear partial_y = f.derivative(1);
    TrivariateLinear partial_z = f.derivative(2);

    // 直交条件 ∇f(x,y,z) * τ
    TrivariateLinear dot;
    dot.coeffs =
    - camera(0) * partial_x.coeffs +
    - camera(1) * partial_y.coeffs +
    - camera(2) * partial_z.coeffs;
    return dot;
}


inline TrivariateLinear perspective_orthogonality_condition(TrivariateQuadratic const& f,
                                                            Eigen::Vector3d const& camera)
{
    Eigen::Matrix3d A;
    Eigen::Vector3d b;
    double c;
    A <<
    f.coeffs(0),     f.coeffs(3)/2.0, f.coeffs(5)/2.0,
    f.coeffs(3)/2.0, f.coeffs(1),      f.coeffs(4)/2.0,
    f.coeffs(5)/2.0, f.coeffs(4)/2.0, f.coeffs(2);
    b << f.coeffs(6)/2.0, f.coeffs(7)/2.0, f.coeffs(8)/2.0;
    c = f.coeffs(9);

    // 直交条件 ∇f(x) * (c-x)
    TrivariateLinear dot;
    dot.coeffs.head(3) = (camera.transpose() * A).transpose() + b;
    dot.coeffs(3) = b.dot(camera) + c;

    return dot;
}

enum class UnivariateLinearTerm
{
    x, y, z, constant
};

class RemoveTerm
{
    BivariateLinear term;
    UnivariateLinearTerm strongest;
};


inline BivariateQuadratic convert_to_algebraic_curve(TrivariateQuadratic const& f,
                                                     Eigen::Vector3d const& camera,
                                                     int const& projection_method,
                                                     BivariateLinear& removed_term,
                                                     int& strongest_coeff)
//                                                     RemoveTerm& strongest)
{
    // 直交条件 ∇f(x,y,z)・τ
    TrivariateLinear dot;
    if (projection_method == 0) dot = orthographic_orthogonality_condition(f, camera);
    else dot = perspective_orthogonality_condition(f, camera);

    // 削除する項
    Eigen::Vector3d linear_coeffs(dot.coeffs(0), dot.coeffs(1), dot.coeffs(2));
    linear_coeffs.cwiseAbs().maxCoeff(&strongest_coeff);
    
    removed_term.coeffs(0) = - dot.coeffs((strongest_coeff + 1)%3) / dot.coeffs(strongest_coeff);
    removed_term.coeffs(1) = - dot.coeffs((strongest_coeff + 2)%3) / dot.coeffs(strongest_coeff);
    removed_term.coeffs(2) = - dot.coeffs(3) / dot.coeffs(strongest_coeff);

    // 輪郭線の代数曲線
//    curve.coeffs(0) = f.coeffs((strongest_coeff + 1)%3);
//    curve.coeffs(1) = f.coeffs((strongest_coeff + 2)%3);
//    curve.coeffs(2) = f.coeffs(4 + (strongest_coeff == 0) - (strongest_coeff == 2));
//    curve.coeffs(3) = f.coeffs(7 + (strongest_coeff == 0) - (strongest_coeff == 2));
//    curve.coeffs(4) = f.coeffs((strongest_coeff + 2)%3 + 6);
//    curve.coeffs(5) = f.coeffs(9);
//
//    curve.coeffs +=
//    f.coeffs(strongest_coeff) * removed_term.square().coeffs +                       // x^2, y^2, z^2
//    f.coeffs((strongest_coeff + 2)%3 + 3) * removed_term.multiple_y().coeffs +       // zx, xy, yz
//    f.coeffs(strongest_coeff + 3) * removed_term.multiple_x().coeffs +               // xy, yz, zx
//    f.coeffs(strongest_coeff + 6) * removed_term.convert_to_quadratic().coeffs;      // x, y, z

    BivariateQuadratic curve;
    if (strongest_coeff == 2) {     // zの項 != 0
        curve.coeffs(0) = f.coeffs(0);
        curve.coeffs(1) = f.coeffs(1);
        curve.coeffs(2) = f.coeffs(3);
        curve.coeffs(3) = f.coeffs(6);
        curve.coeffs(4) = f.coeffs(7);
        curve.coeffs(5) = f.coeffs(9);

        curve.coeffs +=
        f.coeffs(2) * removed_term.square().coeffs +                       // f_2 * z^2
        f.coeffs(4) * removed_term.multiple_y().coeffs +                   // f_4 * yz
        f.coeffs(5) * removed_term.multiple_x().coeffs +                   // f_5 * zx
        f.coeffs(8) * removed_term.convert_to_quadratic().coeffs;          // f_8 * z
    } else if (strongest_coeff == 1) {  // yの項 != 0
        curve.coeffs(0) = f.coeffs(2);
        curve.coeffs(1) = f.coeffs(0);
        curve.coeffs(2) = f.coeffs(5);
        curve.coeffs(3) = f.coeffs(8);
        curve.coeffs(4) = f.coeffs(6);
        curve.coeffs(5) = f.coeffs(9);

        curve.coeffs +=
        f.coeffs(1) * removed_term.square().coeffs +                       // f_1 * y^2
        f.coeffs(3) * removed_term.multiple_y().coeffs +                   // f_3 * xy
        f.coeffs(4) * removed_term.multiple_x().coeffs +                   // f_4 * yz
        f.coeffs(7) * removed_term.convert_to_quadratic().coeffs;          // f_7 * y
    } else if (strongest_coeff == 0) {  // xの項 != 0
        curve.coeffs(0) = f.coeffs(1);
        curve.coeffs(1) = f.coeffs(2);
        curve.coeffs(2) = f.coeffs(4);
        curve.coeffs(3) = f.coeffs(7);
        curve.coeffs(4) = f.coeffs(8);
        curve.coeffs(5) = f.coeffs(9);

        curve.coeffs +=
        f.coeffs(0) * removed_term.square().coeffs +                       // f_0 * x^2
        f.coeffs(5) * removed_term.multiple_y().coeffs +                   // f_3 * zx
        f.coeffs(3) * removed_term.multiple_x().coeffs +                   // f_5 * xy
        f.coeffs(6) * removed_term.convert_to_quadratic().coeffs;          // f_6 * x
    }

    if (coeff_equal_zero(2, dot.coeffs) &&
        coeff_equal_zero(1, dot.coeffs) &&
        coeff_equal_zero(0, dot.coeffs)) {
        strongest_coeff = -1;
    }
    return curve;
}



inline void compute_contour(Eigen::MatrixXd const& P,
                            Eigen::MatrixXi const& Tet,
                            std::vector<TrivariateQuadratic> const& implicit_surfaces,
                            Eigen::Vector3d const& camera,
                            int const& projection_method,
                            std::vector<Eigen::Vector<double, 10>> const& bezier_coeffs,
                            Curve& curve)
{
    int non_rejected_tet_num = 0;
    for (int i = 0; i < implicit_surfaces.size(); i++) {

        if (is_trivial_reject_tet(bezier_coeffs[i])) continue;
        non_rejected_tet_num ++;
        // 四面体の頂点・法線
        std::vector<Eigen::Vector3d> point(4);
        std::vector<Eigen::Vector3d> normal(4);
        tetrahedral_faces(P, Tet, i, point, normal);

        // 代数曲線の計算
        BivariateQuadratic algebraic_curve;
        BivariateLinear removed_term;
        int strongest_coeff;     // z : 2, y : 1, x : 0, invalid : -1
        algebraic_curve = convert_to_algebraic_curve(implicit_surfaces[i], camera, projection_method, removed_term, strongest_coeff);
        if (strongest_coeff == -1) continue;
        // パラメータ化
        ParametrizeConic pc;
        pc.identify_conic(algebraic_curve);     // 二次曲線の識別
        if (pc.conic_type == ConicType::invalid) continue;
        // トリミング
        auto interval = trimming_conic(point, normal, pc, removed_term, strongest_coeff);
        // サンプリング
        pc.parametrize_conic(removed_term, strongest_coeff, interval, curve);
    }
    std::cout << "non rejected tet " << non_rejected_tet_num << std::endl;
    std::cout << "non rejected tet rate " << (double)non_rejected_tet_num/(double)Tet.rows() << std::endl;
}




void test(Eigen::MatrixXd const& P,
          Eigen::MatrixXi const& Tet,
          TrivariateQuadratic const& implicit_surface,
          Eigen::Vector3d const& camera,
          int const& projection_method,
          std::vector<Eigen::Vector<double, 10>> const& bezier_coeffs,
          Curve& curve)
{
//    int non_rejected_tet_num = 0;
    for (int i = 0; i < Tet.rows(); i++)
    {
//        if (is_trivial_reject_tet(bezier_coeffs[i])) continue;
//        non_rejected_tet_num ++;
        // 四面体の頂点・法線
        std::vector<Eigen::Vector3d> point(4);
        std::vector<Eigen::Vector3d> normal(4);
        tetrahedral_faces(P, Tet, i, point, normal);

        // 代数曲線の計算
        BivariateQuadratic algebraic_curve;
        BivariateLinear removed_term;
        int strongest_coeff;     // z : 2, y : 1, x : 0, invalid : -1
        algebraic_curve = convert_to_algebraic_curve(implicit_surface, camera, projection_method, removed_term, strongest_coeff);
        if (strongest_coeff == -1) continue;
        // パラメータ化
        ParametrizeConic pc;
        std::string conic_type = "invalid";
        // 二次曲線の識別
        pc.identify_conic(algebraic_curve);
        // トリミング
        auto interval = trimming_conic(point, normal, pc, removed_term, strongest_coeff);
        // サンプリング
        pc.parametrize_conic(removed_term, strongest_coeff, interval, curve);
    }

//    std::cout << "non rejected tet " << non_rejected_tet_num << std::endl;
//    std::cout << "non rejected tet rate " << (double)non_rejected_tet_num/(double)Tet.rows() << std::endl;
}

#endif /* compute_contour_h */
