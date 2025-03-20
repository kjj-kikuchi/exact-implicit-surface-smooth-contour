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

// 平行投影
inline TrivariateLinear orthographic_contour(TrivariateQuadratic const& f,
                                             Eigen::Vector3d const& camera)
{
    // 勾配ベクトル ∇f(x,y,z)
    TrivariateLinear partial_x = f.derivative(0);
    TrivariateLinear partial_y = f.derivative(1);
    TrivariateLinear partial_z = f.derivative(2);

    // 直交条件 ∇f(x,y,z)・τ
    TrivariateLinear dot;
    dot.coeffs =
    - camera(0) * partial_x.coeffs +
    - camera(1) * partial_y.coeffs +
    - camera(2) * partial_z.coeffs;
    return dot;
}

// 透視投影
inline TrivariateLinear perspective_contour(TrivariateQuadratic const& f,
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

    // 直交条件 ∇f(x)・(c-x)
    TrivariateLinear dot;
    dot.coeffs.head(3) = (camera.transpose() * A).transpose() + b;
    dot.coeffs(3) = b.dot(camera) + c;

    return dot;
}


inline BivariateQuadratic compute_algebraic_curve(TrivariateQuadratic const& f,
                                                  Eigen::Vector3d const& camera,
                                                  int const& projection_method,
                                                  BivariateLinear& removed_term,
                                                  int& strongest_coeff)
{
    // 直交条件 ∇f(x,y,z)・τ
    TrivariateLinear dot;
    if (projection_method == 0) dot = orthographic_contour(f, camera);
    else dot = perspective_contour(f, camera);

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
    if (strongest_coeff == 2) // zの項 != 0
    {
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
    } else if (strongest_coeff == 1) // yの項 != 0
    {
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
    } else if (strongest_coeff == 0) // xの項 != 0
    {
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


//    std::cout << "implicit " << f.coeffs.transpose() << std::endl;
//    std::cout << "partial_x " << partial_x.coeffs.transpose() << std::endl;
//    std::cout << "partial_y " << partial_y.coeffs.transpose() << std::endl;
//    std::cout << "partial_z " << partial_z.coeffs.transpose() << std::endl;
//    std::cout << "dot " << dot.coeffs.transpose() << std::endl;
//    std::cout << "removed_term " << removed_term.coeffs.transpose() << std::endl;
//    std::cout << "curve " << curve.coeffs.transpose() << std::endl;
//    if (coeff_equal_zero(2, dot.coeffs)) std::cout << "!\n";
    if (coeff_equal_zero(2, dot.coeffs) && coeff_equal_zero(1, dot.coeffs) && coeff_equal_zero(0, dot.coeffs))
    {
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
    Eigen::Vector<int, 7> conic_type_num;
    conic_type_num.setZero();
    int non_rejected_tet_num = 0;
    for (int i = 0; i < implicit_surfaces.size(); i++)
    {
        if (is_trivial_reject_tet(bezier_coeffs[i])) continue;
        non_rejected_tet_num ++;
        // 四面体の頂点・法線
        std::vector<Eigen::Vector3d> point(4);
        std::vector<Eigen::Vector3d> normal(4);
        for (int j = 0; j < 4; j++) point[j] = P.row(Tet(i, j));
        normal[0] = get_normal(P.row(Tet(i, 0)), P.row(Tet(i, 1)), P.row(Tet(i, 3)));
        normal[1] = get_normal(P.row(Tet(i, 0)), P.row(Tet(i, 2)), P.row(Tet(i, 1)));
        normal[2] = get_normal(P.row(Tet(i, 0)), P.row(Tet(i, 3)), P.row(Tet(i, 2)));
        normal[3] = get_normal(P.row(Tet(i, 1)), P.row(Tet(i, 2)), P.row(Tet(i, 3)));

        // 代数曲線の計算
        BivariateQuadratic algebraic_curve;
        BivariateLinear removed_term;
        int strongest_coeff;     // z : 2, y : 1, x : 0, invalid : -1
        algebraic_curve = compute_algebraic_curve(implicit_surfaces[i], camera, projection_method, removed_term, strongest_coeff);
        if (strongest_coeff == -1) continue;
        // パラメータ化
        ParametrizeConic pc;
        std::string conic_type = "invalid";
        pc.identify_conic(algebraic_curve, conic_type);     // 二次曲線の識別
//        std::cout << conic_type << std::endl;
        if (conic_type == "ellipse") conic_type_num(0) ++;
        if (conic_type == "hyperbola_1") conic_type_num(1) ++;
        if (conic_type == "hyperbola_2") conic_type_num(2) ++;
        if (conic_type == "intersecting_lines") conic_type_num(3) ++;
        if (conic_type == "parabola") conic_type_num(4) ++;
        if (conic_type == "parallel_lines") conic_type_num(5) ++;
        if (conic_type == "invalid"){
            conic_type_num(6) ++;
            continue;
        }
        // トリミング
        std::array<std::vector<Eigen::Vector2d>, 2> interval;
        interval = trimming_conic(point, normal, pc, removed_term, strongest_coeff, conic_type);
        // サンプリング
        pc.parametrize_conic(removed_term, strongest_coeff, interval, conic_type, curve);
    }
    std::cout << std::endl;
    std::cout << "ellipse " << conic_type_num(0) << std::endl;
    std::cout << "hyperbola_1 " << conic_type_num(1) << std::endl;
    std::cout << "hyperbola_2 " << conic_type_num(2) << std::endl;
    std::cout << "intersecting_lines " << conic_type_num(3) << std::endl;
    std::cout << "parabola " << conic_type_num(4) << std::endl;
    std::cout << "parallel_lines " << conic_type_num(5) << std::endl;
    std::cout << "invalid " << conic_type_num(6) << std::endl;
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
        for (int j = 0; j < 4; j++) point[j] = P.row(Tet(i, j));
        normal[0] = get_normal(P.row(Tet(i, 0)), P.row(Tet(i, 1)), P.row(Tet(i, 3)));
        normal[1] = get_normal(P.row(Tet(i, 0)), P.row(Tet(i, 2)), P.row(Tet(i, 1)));
        normal[2] = get_normal(P.row(Tet(i, 0)), P.row(Tet(i, 3)), P.row(Tet(i, 2)));
        normal[3] = get_normal(P.row(Tet(i, 1)), P.row(Tet(i, 2)), P.row(Tet(i, 3)));

        // 代数曲線の計算
        BivariateQuadratic algebraic_curve;
        BivariateLinear removed_term;
        int strongest_coeff;     // z : 2, y : 1, x : 0, invalid : -1
        algebraic_curve = compute_algebraic_curve(implicit_surface, camera, projection_method, removed_term, strongest_coeff);
        if (strongest_coeff == -1) continue;
        // パラメータ化
        ParametrizeConic pc;
        std::string conic_type = "invalid";
        // 二次曲線の識別
        pc.identify_conic(algebraic_curve, conic_type);
        // トリミング
        auto interval = trimming_conic(point, normal, pc, removed_term, strongest_coeff, conic_type);
        // サンプリング
        pc.parametrize_conic(removed_term, strongest_coeff, interval, conic_type, curve);
    }

//    std::cout << "non rejected tet " << non_rejected_tet_num << std::endl;
//    std::cout << "non rejected tet rate " << (double)non_rejected_tet_num/(double)Tet.rows() << std::endl;
}

#endif /* compute_contour_h */



////一つの四面体と入力メッシュに対して計算
//void test(Eigen::Matrix<double, 10, 3>& P,
//          Eigen::Vector<int, 10>& control,
//          Eigen::Vector<double, 10>& S,
//          Eigen::MatrixXd const& V,
//          Eigen::MatrixXi const& F,
//          Eigen::Vector3d const& camera,
//          int const& projection_method,
//          Curve& curve)
//{
//    P.row(0) = Eigen::Vector3d{ 0.3,  0.3,    0};
//    P.row(1) = Eigen::Vector3d{ 0.3, -0.3,  0.3};
//    P.row(2) = Eigen::Vector3d{ 0.3, -0.3, -0.3};
//    P.row(3) = Eigen::Vector3d{-0.3, -0.3,    0};
//
//    P.row(4) = (P.row(0) + P.row(1))/2.0;
//    P.row(5) = (P.row(0) + P.row(2))/2.0;
//    P.row(6) = (P.row(0) + P.row(3))/2.0;
//    P.row(7) = (P.row(1) + P.row(2))/2.0;
//    P.row(8) = (P.row(1) + P.row(3))/2.0;
//    P.row(9) = (P.row(2) + P.row(3))/2.0;
//
//    control << 0, 1, 2, 3, 4, 5, 6, 7, 8, 9;
//
//    Eigen::VectorXi I;
//    Eigen::MatrixXd C;
//    Eigen::MatrixXd N;
//    igl::signed_distance(P, V, F, igl::SIGNED_DISTANCE_TYPE_PSEUDONORMAL, S, I, C, N);
//
//
//    TrivariateQuadratic implicit_surface;
//
//    Eigen::Matrix<double, 10, 10> Mat;
//    for (int i = 0; i < 10; i++)
//    {
//        Eigen::Vector3d x = P.row(control(i));
//        Mat.row(i) <<
//        x(0)*x(0), x(1)*x(1), x(2)*x(2),
//        x(0)*x(1), x(1)*x(2), x(2)*x(0),
//        x(0), x(1), x(2), 1;
//    }
//    implicit_surface.coeffs = Mat.inverse() * S;
//
//    BivariateQuadratic algebraic_curve;
//    BivariateLinear removed_term;
//    int strongest_coeff;
//    algebraic_curve = compute_algebraic_curve(implicit_surface, camera, projection_method, removed_term, strongest_coeff);
//
//    std::vector<Eigen::Vector3d> point(4);
//    std::vector<Eigen::Vector3d> normal(4);
//    for (int j = 0; j < 4; j++) point[j] = P.row(j);
//    normal[0] = get_normal(P.row(0), P.row(1), P.row(2));
//    normal[1] = get_normal(P.row(0), P.row(3), P.row(1));
//    normal[2] = get_normal(P.row(0), P.row(2), P.row(3));
//    normal[3] = get_normal(P.row(1), P.row(2), P.row(3));
//
//
//    ParametrizeConic pc;
//    std::string conic_type;
//    pc.identify_conic(algebraic_curve, conic_type);
//
//    std::array<std::vector<Eigen::Vector2d>, 2> interval = trimming_conic(point, normal, pc, removed_term, strongest_coeff, conic_type);
//    std::cout << conic_type << std::endl;
//    std::cout << "result\n";
//    for (int i = 0; i < 2; i++) {
//        for (int j = 0; j < interval[i].size(); j++) {
//            std::cout << i << j <<  " [" << interval[i][j](0) << ", " << interval[i][j](1) << "]" << std::endl;
//        }
//    }
//    if (! conic_type.empty()) pc.parametrize_conic(removed_term, strongest_coeff, interval, conic_type, curve);
//}



//if (i == 474 || i == 475)
//        {
//            std::cout << "i " << i << std::endl;
//            std::cout << "implicit " << implicit_surfaces[i].coeffs.transpose() << std::endl;
//            std::cout << "removed term " << removed_term.coeffs.transpose() << std::endl;
//            std::cout << "strongest coeff " << strongest_coeff << std::endl;
//            std::cout << "curve " << algebraic_curve.coeffs.transpose() << std::endl;
//            std::cout << "A\n" << pc.A << std::endl;
//            std::cout << "b " << pc.b.transpose() << std::endl;
//            std::cout << "c " << pc.c << std::endl;
//            std::cout << "U\n" << pc.U << std::endl;
//            std::cout << "sigma " << pc.sigma.transpose() << std::endl;
//            std::cout << "conic type " << conic_type << std::endl;
//            std::cout << "c_hat " << pc.c_hat << std::endl;
//            std::cout << "k " << pc.k.transpose() << std::endl;
//            std::cout << "b_hat " << pc.b_hat.transpose() << std::endl;
//        }
