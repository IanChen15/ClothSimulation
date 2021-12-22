#include <Eigen/Dense>
#include <EigenTypes.h>

#include <algorithm>
void compute_bend_condition(double &bend_condition, 
                            Eigen::Ref<const Eigen::Vector9d> traingle_verticles_1, 
                            Eigen::Ref<const Eigen::Vector9d> traingle_verticles_2);

void d_triangle_norm_dx(Eigen::Matrix<double, 3, 12>  &result, 
                        Eigen::Ref<const Eigen::Vector9d> triangle_vertices);

void d_triangle2_norm_dx(Eigen::Matrix<double, 3, 12>  &result, 
                         Eigen::Ref<const Eigen::Vector9d> triangle_vertices);

void compute_d_bend_condition_dx(Eigen::Vector12d &d_bend_condition_dx,
                                 Eigen::Ref<const Eigen::Vector9d> triangle_vertices_1, 
                                 Eigen::Ref<const Eigen::Vector9d> triangle_vertices_2);

void d2_triangle_norm_dx2(Eigen::Matrix<double, 36, 12> &result, 
                          Eigen::Ref<const Eigen::Vector9d> triangle_vertices);

void d2_triangle2_norm_dx2(Eigen::Matrix<double, 36, 12> &result, 
                           Eigen::Ref<const Eigen::Vector9d> triangle_vertices);
void compute_d2_bend_condition_dx2(Eigen::Matrix1212d &d2_bend_condition_dx2,
                                   Eigen::Ref<const Eigen::Vector9d> triangle_vertices_1, 
                                   Eigen::Ref<const Eigen::Vector9d> triangle_vertices_2);
