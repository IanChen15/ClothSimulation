#include <compute_stretch_condition.h>

#include<iostream>

void stretch_derivative(Eigen::Matrix39d &stretch_derivative,
                        double u, double v);

void compute_stretch(Eigen::Matrix32d &stretch, 
                     Eigen::Ref<const Eigen::Vector6d> ref_plane_coordinate, 
                     Eigen::Ref<const Eigen::Vector9d> triangle_vertices) {

    Eigen::Vector3d xi = triangle_vertices.segment<3>(0);
    Eigen::Vector3d xj = triangle_vertices.segment<3>(3);
    Eigen::Vector3d xk = triangle_vertices.segment<3>(6);

    Eigen::Vector2d pi = ref_plane_coordinate.segment<2>(0);
    Eigen::Vector2d pj = ref_plane_coordinate.segment<2>(2);
    Eigen::Vector2d pk = ref_plane_coordinate.segment<2>(4);
    
    Eigen::Matrix32d dx_matrix;
    dx_matrix.col(0) = xj - xi;
    dx_matrix.col(1) = xk - xi;

    Eigen::Matrix2d dp_matrix;
    dp_matrix.col(0) = pj - pi;
    dp_matrix.col(1) = pk - pi;

    stretch = dx_matrix * dp_matrix.inverse();
};

void compute_stretch_condition(Eigen::Vector2d &stretch_condition, 
                               Eigen::Ref<const Eigen::Matrix32d> stretch, 
                               double area, double bu, double bv) {
    stretch_condition(0) = area * (stretch.col(0).norm() - bu);
    stretch_condition(1) = area * (stretch.col(1).norm() - bv);
};

void d_stretch_dx(Eigen::Matrix39d &d_stretch_u_dx, 
                  Eigen::Matrix39d &d_stretch_v_dx,
                  Eigen::Ref<const Eigen::Vector6d> ref_plane_coordinate) {
    Eigen::Vector2d pi = ref_plane_coordinate.segment<2>(0);
    Eigen::Vector2d pj = ref_plane_coordinate.segment<2>(2);
    Eigen::Vector2d pk = ref_plane_coordinate.segment<2>(4);
    
    Eigen::Matrix2d dp_matrix;
    dp_matrix.col(0) = pj - pi;
    dp_matrix.col(1) = pk - pi;
    Eigen::Matrix2d dp_inv_matrix = dp_matrix.inverse();
    
    stretch_derivative(d_stretch_u_dx, dp_inv_matrix(0, 0), dp_inv_matrix(1, 0));
    stretch_derivative(d_stretch_v_dx, dp_inv_matrix(0, 1), dp_inv_matrix(1, 1));    
}

void stretch_derivative(Eigen::Matrix39d &stretch_derivative,
                        double u, double v) {
    stretch_derivative << -u-v,    0,    0, u, 0, 0, v, 0, 0, 
                             0, -u-v,    0, 0, u, 0, 0, v, 0,
                             0,    0, -u-v, 0, 0, u, 0, 0, v;
}


void compute_d_stretch_condition_dx(Eigen::Matrix<double, 2, 9> &d_stretch_condition_dx, 
                                    Eigen::Ref<const Eigen::Matrix39d> d_stretch_u_dx, 
                                    Eigen::Ref<const Eigen::Matrix39d> d_stretch_v_dx,
                                    Eigen::Ref<const Eigen::Matrix32d> stretch, 
                                    double area) {

    d_stretch_condition_dx.row(0) = 
        area * d_stretch_u_dx.transpose() * stretch.col(0) / stretch.col(0).norm();

    d_stretch_condition_dx.row(1) = 
        area * d_stretch_v_dx.transpose() * stretch.col(1) / stretch.col(1).norm();
}

void compute_d2_stretch_condition_dx2(Eigen::Matrix99d &d2_stretch_condition_u_dx2, 
                                      Eigen::Matrix99d &d2_stretch_condition_v_dx2, 
                                      Eigen::Ref<const Eigen::Matrix39d> d_stretch_u_dx, 
                                      Eigen::Ref<const Eigen::Matrix39d> d_stretch_v_dx,
                                      Eigen::Ref<const Eigen::Matrix32d> stretch, 
                                      double area) {
    Eigen::Vector3d wu = stretch.col(0); 
    Eigen::Vector3d wv = stretch.col(1); 

    double wu_norm = wu.norm();
    double wv_norm = wv.norm();
    Eigen::Matrix93d temp_wu = d_stretch_u_dx.transpose() * (Eigen::Matrix3d::Identity() - wu * wu.transpose() / (wu_norm * wu_norm)) / wu_norm;
    Eigen::Matrix93d temp_wv = d_stretch_v_dx.transpose() * (Eigen::Matrix3d::Identity() - wv * wv.transpose() / (wv_norm * wv_norm)) / wv_norm;
    
    d2_stretch_condition_u_dx2 = area * temp_wu * d_stretch_u_dx;
    d2_stretch_condition_v_dx2 = area * temp_wv * d_stretch_v_dx;
}