#include <Eigen/Dense>
#include <EigenTypes.h>

void compute_shear_condition(double &shear_condition, 
                               Eigen::Ref<const Eigen::Matrix32d> stretch_force, 
                               double area);

void compute_d_shear_condition_dx(Eigen::Vector9d &d_shear_condition_dx, 
                                  Eigen::Ref<const Eigen::Matrix39d> d_stretch_u_dx, 
                                  Eigen::Ref<const Eigen::Matrix39d> d_stretch_v_dx,
                                  Eigen::Ref<const Eigen::Matrix32d> stretch_force, 
                                  double area);

void compute_d2_shear_condition_dx2(Eigen::Matrix99d &d2_shear_condition_dx2, 
                                    Eigen::Ref<const Eigen::Matrix39d> d_stretch_u_dx, 
                                    Eigen::Ref<const Eigen::Matrix39d> d_stretch_v_dx,
                                    double area);
