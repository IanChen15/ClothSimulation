#include <Eigen/Dense>
#include <EigenTypes.h>

void compute_stretch(Eigen::Matrix32d &stretch, 
                     Eigen::Ref<const Eigen::Vector6d> ref_plane_coordinate, 
                     Eigen::Ref<const Eigen::Vector9d> traingle_verticles);

void compute_stretch_condition(Eigen::Vector2d &stretch_condition, 
                               Eigen::Ref<const Eigen::Matrix32d> stretch, 
                               double area, double bu, double bv);

void d_stretch_dx(Eigen::Matrix39d &d_stretch_u_dx, 
                  Eigen::Matrix39d &d_stretch_v_dx,
                  Eigen::Ref<const Eigen::Vector6d> ref_plane_coordinate);

void compute_d_stretch_condition_dx(Eigen::Matrix<double, 2, 9> &d_stretch_condition_dx, 
                                    Eigen::Ref<const Eigen::Matrix39d> d_stretch_u_dx, 
                                    Eigen::Ref<const Eigen::Matrix39d> d_stretch_v_dx,
                                    Eigen::Ref<const Eigen::Matrix32d> stretch, 
                                    double area);
void compute_d2_stretch_condition_dx2(Eigen::Matrix99d &d2_stretch_condition_u_dx2, 
                                      Eigen::Matrix99d &d2_stretch_condition_v_dx2, 
                                      Eigen::Ref<const Eigen::Matrix39d> d_stretch_u_dx, 
                                      Eigen::Ref<const Eigen::Matrix39d> d_stretch_v_dx,
                                      Eigen::Ref<const Eigen::Matrix32d> stretch, 
                                      double area);