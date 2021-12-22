#include <compute_shear_condition.h>
#include<iostream>

void compute_shear_condition(double &shear_condition, 
                             Eigen::Ref<const Eigen::Matrix32d> stretch_force, 
                             double area) {

    shear_condition = area * stretch_force.col(0).dot(stretch_force.col(1));
};


void compute_d_shear_condition_dx(Eigen::Vector9d &d_shear_condition_dx, 
                                  Eigen::Ref<const Eigen::Matrix39d> d_stretch_u_dx, 
                                  Eigen::Ref<const Eigen::Matrix39d> d_stretch_v_dx,
                                  Eigen::Ref<const Eigen::Matrix32d> stretch_force, 
                                  double area) {
    d_shear_condition_dx = area * (d_stretch_v_dx.transpose() * stretch_force.col(0)  
                                    + d_stretch_u_dx.transpose() * stretch_force.col(1));
}

void compute_d2_shear_condition_dx2(Eigen::Matrix99d &d2_shear_condition_dx2, 
                                    Eigen::Ref<const Eigen::Matrix39d> d_stretch_u_dx, 
                                    Eigen::Ref<const Eigen::Matrix39d> d_stretch_v_dx,
                                    double area) {
    d2_shear_condition_dx2 = area * (d_stretch_v_dx.transpose() * d_stretch_u_dx
                                    + d_stretch_u_dx.transpose() * d_stretch_v_dx);
}