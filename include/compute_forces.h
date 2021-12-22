#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <EigenTypes.h>

void compute_stretch_force(Eigen::Vector9d &stretch_force, 
                           Eigen::Matrix99d &d_stretch_force_dx,
                           Eigen::Vector9d &stretch_damping_force, 
                           Eigen::Matrix99d &d_stretch_damping_force_dx,
                           Eigen::Matrix99d &d_stretch_damping_force_dv,
                           Eigen::Ref<const Eigen::Vector9d> triangle_dot,
                           Eigen::Ref<const Eigen::Matrix32d> stretch,
                           Eigen::Ref<const Eigen::Matrix39d> d_stretch_u_dx,
                           Eigen::Ref<const Eigen::Matrix39d> d_stretch_v_dx,
                           double k, double k_damp,
                           double area, double bu, double bv);
void compute_shear_force(Eigen::Vector9d &shear_force, 
                         Eigen::Matrix99d &d_shear_force_dx, 
                         Eigen::Vector9d &shear_damping_force,
                         Eigen::Matrix99d &d_shear_damping_force_dx,
                         Eigen::Matrix99d &d_shear_damping_force_dv,
                         Eigen::Ref<Eigen::Vector9d> triangle_dot,
                         Eigen::Ref<const Eigen::Matrix32d> stretch,
                         Eigen::Ref<const Eigen::Matrix39d> d_stretch_u_dx,
                         Eigen::Ref<const Eigen::Matrix39d> d_stretch_v_dx,
                         double k, double k_damp, double area);
void compute_bend_force(Eigen::Vector12d &bend_force,
                        Eigen::Matrix1212d &d_bend_force_dx, 
                        Eigen::Vector12d &bend_damping_force,
                        Eigen::Matrix1212d &d_bend_damping_force_dx,
                        Eigen::Matrix1212d &d_bend_damping_force_dv,
                        Eigen::Ref<const Eigen::Vector9d> triangle_dot1, 
                        Eigen::Ref<const Eigen::Vector9d> triangle_dot2, 
                        Eigen::Ref<const Eigen::Vector9d> ordered_triangle1, 
                        Eigen::Ref<const Eigen::Vector9d> ordered_triangle2, 
                        double k, double k_damp);