//
//  main.cpp
//  piecewise-quadraic-implicit-surface
//
//  Created by 菊池祐作 on 2024/12/23.
//

#include <igl/readOBJ.h>
#include "bcc_lattice.h"
#include "compute_contour.h"
#include "export_file.h"


int main(int argc, const char * argv[])
{
    if (argc != 2)
    {
        std::cout << "wrong command line argument" << std::endl;
        std::exit(1);
    }
    Eigen::Vector3d camera;
    camera << 0, 0.5, 1.5;
//    camera << 2, 2, 2;
//    camera << 1, -2, 0;
//    camera << -1, -1, -1;
//    camera << -1, 1, 1;
    double cell_resolution = 0.02;
    int projection_method = 1;  // orthographic : 0, perspective : 1


    // Read Mesh...........................................................
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    std::string filename = std::string(argv[1]);
    igl::readOBJ(filename, V, F);

    std::string meshname = filename;
    meshname = meshname.erase(meshname.length()-4) + "_";

//    auto start = std::chrono::system_clock::now();


    // Normalize mesh.........................................................
    Eigen::Vector3d min = V.colwise().minCoeff();
    Eigen::Vector3d max = V.colwise().maxCoeff();
    double scale = 1.0 / (max - min).norm();
    Eigen::Vector3d center = (max + min) * 0.5;
    for (int i = 0; i < V.rows(); i++) V.row(i) = (V.row(i) - center.transpose()) * scale;


    // Construct BCC lattice.................................................
    Eigen::MatrixXd P;                                      // 格子点
    Eigen::MatrixXd midP;                                   // 辺の中点
    std::map<std::pair<int, int>, int> midP_idx;            // 辺のペアから辺の中点の添え字へのマップ
    Eigen::MatrixXi Tet;                                    // 四面体リスト
    std::vector<Eigen::Vector<int, 10>> control_idx;        // 四面体の制御点のインデックス
    Eigen::Vector3i n;                                      // 格子点数

    construct_lattice(V, cell_resolution, n, P, Tet);
    enumetate_mid_points(P, Tet, midP, midP_idx, control_idx);
    std::cout << "n " << n.transpose() << std::endl;

    
    // Compute piecewise quadric implicit surface.....................................
    Eigen::VectorXd GS;                                     // BCCのsigned distance
    Eigen::VectorXd midS;                                   // 辺中点のsigned diatance
    std::vector<Eigen::Vector<double, 10>> bezier_coeffs;   // 制御点のベジェ係数
    std::vector<TrivariateQuadratic> implicit_surfaces;     // 四面体の関数f(x,y,z)

    get_signed_distance(P, V, F, GS);
    get_signed_distance(midP, V, F, midS);
    Eigen::VectorXd S(GS.size() + midS.size());
    S << GS, midS;
    compute_implicit_surface(P, S, Tet, control_idx, implicit_surfaces);
    compute_bezier_coefficients(Tet, GS, midS, midP_idx, bezier_coeffs);


    // Compute contour........................................................
    auto start_contour = std::chrono::system_clock::now();
    Curve curve;
    compute_contour(P, Tet, implicit_surfaces, camera, projection_method, bezier_coeffs, curve);
    auto end_contour = std::chrono::system_clock::now();


    // TEST
//    Curve curve2;
//    TrivariateQuadratic implicit_surface;
//    double r = 0.288675;
//    implicit_surface.coeffs = {1, 2, 1, 0, 0, 0, 0, 0, 0, -r*r};    // ellipse
//    implicit_surface.coeffs = {1, 1, 1, 0, 0, 0, 0, 0, 0, -r*r};    // ellipse
//    implicit_surface.coeffs = {1, -1, 1, 0, 0, 0, 0, 0, 0, -0.01};   // hyperbola
//    implicit_surface.coeffs = {1, -2, 1, 0, 0, 0, 0, 0, 0, 0};      // intersecting line
//    implicit_surface.coeffs = {3, 0, 8, 0, 0, 0, 0, -1, 0, -0.25};      // parabola
//    implicit_surface.coeffs = {1, 0, 1, 0, 0, 0, 0, 0, 0, -0.04};       // parallel lines
//    std::vector<Eigen::Vector<double, 10>> bezier_coeffs;
//    compute_bezier_coefficients(P, midP, Tet, control_idx, implicit_surface, midP_idx, bezier_coeffs);

//    auto start_contour = std::chrono::system_clock::now();
//    test(P, Tet, implicit_surface, camera, projection_method, bezier_coeffs, curve2);
//
//    auto end_contour = std::chrono::system_clock::now();


//    auto end = std::chrono::system_clock::now();





    // Show result............................................................

//    export_file(P, GS, meshname + "lattice");
//    export_obj(V, F, meshname + "normalized");
//    export_curve_obj(curve2, meshname + "contour");
//    export_tet_file(Tet);
//    {
//        Eigen::Vector3i n_;
//        Eigen::Vector3d cl = Eigen::Vector3d::Ones() * 0.01;
//        Eigen::Vector3d min = V.colwise().minCoeff();
//        Eigen::Vector3d max = V.colwise().maxCoeff();
//        cl = cl.cwiseProduct(max - min);
//        n_(0) = std::ceil((max(0) - min(0)) / cl(0)) + 1;
//        n_(1) = std::ceil((max(1) - min(1)) / cl(1)) + 1;
//        n_(2) = std::ceil((max(2) - min(2)) / cl(2)) + 1;
//
//        Eigen::MatrixXd sampleP;
//        Eigen::VectorXd SDF;
//        sampleP.resize(n_(0) * n_(1) * n_(2), 3);
//        int point_idx = 0;
//        for (int iz = 0; iz < n_(2); iz++)
//        {
//            for (int iy = 0; iy < n_(1); iy++)
//            {
//                for (int ix = 0; ix < n_(0); ix++)
//                {
//                    sampleP.row(point_idx) = min + Eigen::Vector3d(ix, iy, iz).cwiseProduct(cl);
//                    point_idx ++;
//                }
//            }
//        }
//        get_signed_distance(sampleP, V, F, SDF);
//        export_file(sampleP, SDF, n_, cl, meshname + "sdf"); 
//    }



    export_povray_files(camera, projection_method, implicit_surfaces, P, Tet, bezier_coeffs, curve, meshname);


    // test
//    export_povray_files(camera, projection_method, implicit_surface, P, Tet, bezier_coeffs, curve2, meshname);

    using namespace std::chrono_literals;
    std::cout << "Tet " << Tet.rows() << std::endl;
    std::cout << "Execute time:\n";
    std::cout << "Compute contour : " << (end_contour - start_contour) / 1.0s << " s\n";
//    std::cout << "Compute : " << (end - start) / 1.0s << " s\n";

    return 0;
}







//    // 主グリッドのSDF
//    Eigen::MatrixXd primeP = P.block(0, 0, n(0)*n(1)*n(2), 3);
//    Eigen::VectorXd primeS = GS.head(n(0)*n(1)*n(2));
//    export_file(primeP, primeS, n, cell_size, meshname + "primeP");
//    // 双対グリッドのSDF
//    Eigen::MatrixXd dualP = P.block(n(0)*n(1)*n(2), 0, n(0)*n(1)*n(2), 3);
//    Eigen::VectorXd dualS = GS.tail(n(0)*n(1)*n(2));
//    export_file(dualP, dualS, n, cell_size, meshname + "dualP");





//    // 陰関数のテスト用
//    Eigen::MatrixXd sampleP;
//    Eigen::VectorXd sampleSDF;
//    Eigen::VectorXd sampleFunc1;
//    Eigen::VectorXd sampleFunc2;
//
//    sampleP.resize(Tet.rows(), 3);
//    for (int i = 0; i < Tet.rows(); i++)
//    {
//        sampleP.row(i) = ( P.row(Tet(i, 0)) + P.row(Tet(i, 1)) + P.row(Tet(i, 2)) + P.row(Tet(i, 3)) ) / 4.0;
//    }
//
//    // f(x,y,z)
//    test(sampleP, coef, sampleFunc1);
//    export_file(sampleP, sampleFunc1, meshname + "func");
//
//    get_signed_distance(sampleP, V, F, sampleSDF);
//    export_file(sampleP, sampleSDF, meshname + "sdf");
//    export_tet_file(Tet);
//    export_file(coef);

//    // f(x,y,z)とSDFの差
//    Eigen::VectorXd diff = (sampleSDF - sampleFunc1).cwiseAbs();
//    export_file(sampleP, diff, meshname + "diff");
//
//    std::cout << "Max diff: " << diff.maxCoeff() << std::endl;
//    std::cout << "Min diff: " << diff.minCoeff() << std::endl;
//    std::cout << "Mean diff: " << diff.mean() << std::endl;
//
//    std::cout << "P: " << P.rows() << std::endl;
//    std::cout << "midP: " << midP.rows() << std::endl;
//    std::cout << "Tet: " << Tet.rows() << std::endl;
