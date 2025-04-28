//
//  export_file.h
//  piecewise-quadraic-implicit-surface
//
//  Created by 菊池祐作 on 2025/03/06.
//

#ifndef export_file_h
#define export_file_h

#include <fstream>
#include <vector>
#include <Eigen/Core>

inline void export_curve_inc(Curve const& curve,
                             std::string const& name)
{
    std::ofstream of;
    std::string filename = name + "curve.inc";
    of.open(filename, std::ios::out);

//    for (int i = 0; i < curve.L.size(); i++)
    for (auto const& line_segment : curve.L)
    {
        of << "sphere_sweep {\n  linear_spline,\n";
        of << "  " << NUM_SEGMENTS + 1 << ",\n";

        for (int j = 0; j < NUM_SEGMENTS; j++)
        {
            Eigen::Vector3d p = curve.P[line_segment[j](0)];
            of << "<" << p(0) << ", " << p(1) << ", " << p(2) << ">, " << "0.004\n";
            if (j == NUM_SEGMENTS - 1)
            {
                p = curve.P[line_segment[j](0) + 1];
                of << "<" << p(0) << ", " << p(1) << ", " << p(2) << ">, " << "0.004\n";
            }
        }
        of << "}\n";
    }
    of.close();
}

inline void export_piecewise_quadratic_inc(std::vector<TrivariateQuadratic> const& implicit_surfaces,
                                           Eigen::MatrixXd const& P,
                                           Eigen::MatrixXi const& Tet,
                                           std::vector<Eigen::Vector<double, 10>> const& bezier_coeffs,
                                           std::string const& name)
{
    std::ofstream of;
    std::string filename = name + "piecewise_quadratic.inc";
    of.open(filename, std::ios::out);
    for (int i = 0; i < Tet.rows(); i++)
    {
        if (is_trivial_reject_tet(bezier_coeffs[i])) continue;
        
        auto is = implicit_surfaces[i].coeffs;
        of << "quadric {\n";
        of << "  <" << is(0) << ", " << is(1) << ", " << is(2) << ">, ";
        of << "  <" << is(3) << ", " << is(5) << ", " << is(4) << ">, ";
        of << "  <" << is(6) << ", " << is(7) << ", " << is(8) << ">, " << is(9) << "\n";
        of << "  TetrahedralDomain(";
        for (int j = 0; j < 4; j++)
        {
            of << "< " << P(Tet(i, j), 0) << ", " << P(Tet(i, j), 1) << ", " << P(Tet(i, j), 2) << ">";
            if (j == 3) of << ")\n}\n";
            else of << ", ";
        }
    }
    of.close();
}

inline void export_piecewise_quadratic_inc(TrivariateQuadratic const& implicit_surface,
                                           Eigen::MatrixXd const& P,
                                           Eigen::MatrixXi const& Tet,
                                           std::vector<Eigen::Vector<double, 10>> const& bezier_coeffs,
                                           std::string const& name)
{
    std::ofstream of;
    std::string filename = name + "piecewise_quadratic.inc";
    of.open(filename, std::ios::out);
    for (int i = 0; i < Tet.rows(); i++)
    {
        if (is_trivial_reject_tet(bezier_coeffs[i])) continue;

        auto is = implicit_surface.coeffs;
        of << "quadric {\n";
        of << "  <" << is(0) << ", " << is(1) << ", " << is(2) << ">, ";
        of << "  <" << is(3) << ", " << is(5) << ", " << is(4) << ">, ";
        of << "  <" << is(6) << ", " << is(7) << ", " << is(8) << ">, " << is(9) << "\n";
        of << "  TetrahedralDomain(";
        for (int j = 0; j < 4; j++)
        {
            of << "< " << P(Tet(i, j), 0) << ", " << P(Tet(i, j), 1) << ", " << P(Tet(i, j), 2) << ">";
            if (j == 3) of << ")\n}\n";
            else of << ", ";
        }
    }
    of.close();
}

inline void export_povray_files(Eigen::Vector3d const& camera,
                                int const& projection_method,
                                std::vector<TrivariateQuadratic> const& implicit_surface,
                                Eigen::MatrixXd const& P,
                                Eigen::MatrixXi const& Tet,
                                std::vector<Eigen::Vector<double, 10>> const& bezier_coeffs,
                                Curve const& curve,
                                std::string const& name)
{
    std::ofstream of;
    std::string filename = name + "scene.pov";
    of.open(filename, std::ios::out);
    // camera
    of << "camera {\n";
    if (projection_method == 0) of << "  orthographic\n";
    else of << "  perspective\n";
    of << "  location <" << camera(0) << ", " << camera(1) << ", " << camera(2) << ">\n";
    of << "  look_at  <0, 0, 0>\n";
    of << "  right -1.0 * x\n  up 0.75 * y\n}\n\n";
    // light source
    of << "light_source {\n  <5, 5, 5>\n  color rgb <1, 1, 1>\n  parallel\n  point_at <0, 0, 0>\n}\n\n";

    of << "#include \"tetrahedral_domain.inc\"\n\n";
    of << "#default {\n  texture {\n    pigment { color rgb <0.9, 0.9, 0.9> }\n    finish { ambient rgb <0.5, 0.5, 0.5> }\n  }\n}\n\n";
//    of << "#default {\n  texture {\n    pigment { color rgb <1, 1, 1> }\n    finish { ambient rgb <1, 1, 1> }\n  }\n}\n\n";
    of << "#include \"" << name << "piecewise_quadratic.inc\"\n\n";
//    of << "#default {\n  texture {\n    pigment { color rgb <1, 0, 0> }\n    finish { ambient rgb <1, 0, 0> }\n  }\n}\n\n"; // 赤
    of << "#default {\n  texture {\n    pigment { color rgb <0, 0, 0> }\n    finish { ambient rgb <0, 0, 0> }\n  }\n}\n\n"; // 黒
    of << "#include \"" << name << "curve.inc\"\n";
    of.close();

    export_piecewise_quadratic_inc(implicit_surface, P, Tet, bezier_coeffs, name);
    export_curve_inc(curve, name);
}

inline void export_povray_files(Eigen::Vector3d const& camera,
                                int const& projection_method,
                                TrivariateQuadratic const& implicit_surfaces,
                                Eigen::MatrixXd const& P,
                                Eigen::MatrixXi const& Tet,
                                std::vector<Eigen::Vector<double, 10>> const& bezier_coeffs,
                                Curve const& curve,
                                std::string const& name)
{
    std::ofstream of;
    std::string filename = name + "scene.pov";
    of.open(filename, std::ios::out);
    // camera
    of << "camera {\n";
    if (projection_method == 0) of << "  orthographic\n";
    else of << "  perspective\n";
    of << "  location <" << camera(0) << ", " << camera(1) << ", " << camera(2) << ">\n";
    of << "  look_at  <0, 0, 0>\n";
    of << "  right -1.0 * x\n  up 0.75 * y\n}\n\n";
    // light source
    of << "light_source {\n  <5, 5, 5>\n  color rgb <1, 1, 1>\n  parallel\n  point_at <0, 0, 0>\n}\n\n";

    of << "#include \"tetrahedral_domain.inc\"\n\n";
    of << "#default {\n  texture {\n    pigment { color rgb <0.9, 0.9, 0.9> }\n    finish { ambient rgb <0.5, 0.5, 0.5> }\n  }\n}\n\n";
    of << "#include \"" << name << "piecewise_quadratic.inc\"\n\n";
    of << "#default {\n  texture {\n    pigment { color rgb <1, 0, 0> }\n    finish { ambient rgb <1, 0, 0> }\n  }\n}\n\n";
    of << "#include \"" << name << "curve.inc\"\n";
    of.close();

    export_piecewise_quadratic_inc(implicit_surfaces, P, Tet, bezier_coeffs, name);
    export_curve_inc(curve, name);
}


inline void export_curve_obj(Curve const& curve,
                             std::string const& name)
{
    std::ofstream of;
    std::string filename = name + ".obj";
    of.open(filename, std::ios::out);
    for (auto const& p : curve.P)
    {
        of << "v " << p.transpose() << std::endl;
    }
    for (auto const& line_segment : curve.L)
    {
        for (auto const& l : line_segment)
        {
            of << "l " << l(0)+1 << " " << l(1)+1 << std::endl;
        }
    }
    of.close();
}

inline void export_obj(Eigen::MatrixXd const& V,
                       Eigen::MatrixXi const& F,
                       std::string const& name)
{
    std::ofstream of;
    std::string filename = name + ".obj";
    of.open(filename, std::ios::out);
    for(int i = 0; i < V.rows(); i++)
    {
        of << "v " << V.row(i) << std::endl;
    }
    for(int i = 0; i < F.rows(); i++)
    {
        of << "f " << F(i,0)+1 << " " << F(i,1)+1 << " " << F(i,2)+1 << std::endl;
    }
    of.close();
}

inline void export_file(Eigen::MatrixXd const& P,
                        Eigen::VectorXd const& S,
                        Eigen::Vector3i const& n,
                        Eigen::Vector3d const& cell_size,
                        std::string name)
{
    std::ofstream of;
    std::string filename = name + ".vtk";
    of.open(filename, std::ios::out);
    of << "# vtk DataFile Version 3.0\n" << "vtk output\n" << "ASCII\n" << "DATASET STRUCTURED_POINTS\n";
    of << "DIMENSIONS " << n.transpose() << std::endl;
    of << "ORIGIN " << P.row(0) << std::endl;
    of << "SPACING " << cell_size.transpose() << std::endl;
    of << "POINT_DATA " << n(0) * n(1) * n(2) << std::endl;
    of << "SCALARS func_value float 1\n";
    of << "LOOKUP_TABLE default\n";
    for (int i = 0; i < S.size(); i++)
    {
        of << S(i) << std::endl;
    }
    of.close();
}

inline void export_file(Eigen::MatrixXd const& P,
                        Eigen::VectorXd const& S,
                        std::string name)
{
    std::ofstream of;
    std::string filename = name + ".vtk";
    of.open(filename, std::ios::out);
    of << "# vtk DataFile Version 3.0\n" << "vtk output\n" << "ASCII\n" << "DATASET POLYDATA\n";
    of << "POINTS " << P.rows() << " float\n";
    for (int i = 0; i < P.rows(); i++)
    {
        of << P.row(i) << std::endl;
    }
    of << "VERTICES " << P.rows() << " " << 2 * P.rows() << std::endl;
    for (int i = 0; i < P.rows(); i++)
    {
        of << "1 " << i << std::endl;
    }
    of << "POINT_DATA " << S.size() << std::endl;
    of << "SCALARS signed_distance float\n";
    of << "LOOKUP_TABLE default\n";
    for (int i = 0; i < S.size(); i++)
    {
        of << S(i) << std::endl;
    }
    of.close();
}

inline void export_tet_file(Eigen::MatrixXi const& tet)
{
    std::ofstream of;
    std::string filename = "./tet.txt";
    of.open(filename, std::ios::out);
    of << tet << std::endl;
    of.close();
}

inline void export_tet(std::vector<Eigen::Vector3d> const& points, std::string name)
{
    std::ofstream of;
    std::string filename = "./tetra_" + name + ".obj";
    of.open(filename, std::ios::out);
    for(auto const& p : points)
    {
        of << "v " << p.transpose() << std::endl;
    }
    of << "f " << 1 << " " << 2 << " " << 4 << std::endl;
    of << "f " << 1 << " " << 3 << " " << 2 << std::endl;
    of << "f " << 1 << " " << 4 << " " << 3 << std::endl;
    of << "f " << 2 << " " << 3 << " " << 4 << std::endl;
    of.close();
}
#endif /* export_file_h */
