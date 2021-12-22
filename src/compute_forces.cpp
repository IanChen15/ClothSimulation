#include <compute_forces.h>

#include <compute_stretch_condition.h>
#include <compute_shear_condition.h>
#include <compute_bend_condition.h>

#include <iostream>
#define TOL (1e-30)
void remove_small_values(Eigen::Vector9d &vector) {
    for(int i = 0; i < 9; i++) {
        if (std::abs(vector(i)) < TOL) {
            vector(i) = 0;
        }
    }
}
void remove_small_values(Eigen::Vector12d &vector) {
    for(int i = 0; i < 12; i++) {
        if (std::abs(vector(i)) < TOL) {
            vector(i) = 0;
        }
    }
}
void remove_small_values(Eigen::Matrix1212d &m) {
    for(int i = 0; i < 12; i++) {
        for (int j = 0; j < 12; j++) {
            if (std::abs(m(i, j)) < TOL) {
                m(i, j) = 0;
            }
        }
    }
}
void remove_small_values(Eigen::Matrix99d &m) {
    for(int i = 0; i < 9; i++) {
        for (int j = 0; j < 9; j++) {
            if (std::abs(m(i, j)) < TOL) {
                m(i, j) = 0;
            }
        }
    }
}

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
                           double area, double bu, double bv){
    //stretch
    Eigen::Vector2d stretch_condition;
    Eigen::Matrix<double, 2, 9> d_stretch_condition_dx;
    Eigen::Matrix99d d2_stretch_condition_u_dx2;
    Eigen::Matrix99d d2_stretch_condition_v_dx2;

    compute_stretch_condition(stretch_condition, stretch, area, bu, bv);
    compute_d_stretch_condition_dx(d_stretch_condition_dx, 
                                    d_stretch_u_dx, d_stretch_v_dx,
                                    stretch, area);
    compute_d2_stretch_condition_dx2(d2_stretch_condition_u_dx2, d2_stretch_condition_v_dx2, 
                                        d_stretch_u_dx, d_stretch_v_dx, stretch, area);
    stretch_force = k * d_stretch_condition_dx.transpose() * stretch_condition;
    // std::cout << d2_stretch_condition_u_dx2;
    // std::cout << d2_stretch_condition_u_dx2;
    d_stretch_force_dx = k * (
                d_stretch_condition_dx.transpose() * d_stretch_condition_dx 
                + d2_stretch_condition_u_dx2 * stretch_condition(0) 
                + d2_stretch_condition_v_dx2 * stretch_condition(1) 
            );

    Eigen::Vector2d cdot = d_stretch_condition_dx * triangle_dot;
    stretch_damping_force = k_damp * d_stretch_condition_dx.transpose() * cdot;
    d_stretch_damping_force_dx = k_damp * (
                d2_stretch_condition_u_dx2 * cdot(0) 
                + d2_stretch_condition_v_dx2 * cdot(1) 
    );
    d_stretch_damping_force_dv = k_damp * d_stretch_condition_dx.transpose() * d_stretch_condition_dx;

    remove_small_values(stretch_force);
    remove_small_values(stretch_damping_force);
    remove_small_values(d_stretch_force_dx);
    remove_small_values(d_stretch_damping_force_dx);
    remove_small_values(d_stretch_damping_force_dv);
}

void compute_shear_force(Eigen::Vector9d &shear_force, 
                         Eigen::Matrix99d &d_shear_force_dx, 
                         Eigen::Vector9d &shear_damping_force,
                         Eigen::Matrix99d &d_shear_damping_force_dx,
                         Eigen::Matrix99d &d_shear_damping_force_dv,
                         Eigen::Ref<Eigen::Vector9d> triangle_dot,
                         Eigen::Ref<const Eigen::Matrix32d> stretch,
                         Eigen::Ref<const Eigen::Matrix39d> d_stretch_u_dx,
                         Eigen::Ref<const Eigen::Matrix39d> d_stretch_v_dx,
                         double k, double k_damp, double area) {
    double shear_condition;
    Eigen::Vector9d d_shear_condition_dx;
    compute_shear_condition(shear_condition, stretch, area);
    compute_d_shear_condition_dx(d_shear_condition_dx, d_stretch_u_dx, d_stretch_v_dx, stretch, area);
    Eigen::Matrix99d d2_shear_condition_dx2;
    compute_d2_shear_condition_dx2(d2_shear_condition_dx2, d_stretch_u_dx, d_stretch_v_dx, area);

    shear_force = k * d_shear_condition_dx * shear_condition;
    d_shear_force_dx = k * (
                d_shear_condition_dx * d_shear_condition_dx.transpose()
                + shear_condition * d2_shear_condition_dx2
            );

    double c_dot = d_shear_condition_dx.transpose() * triangle_dot;
    shear_damping_force = k_damp * d_shear_condition_dx * c_dot;
    d_shear_damping_force_dx = k_damp * d2_shear_condition_dx2 * c_dot;
    d_shear_damping_force_dv = k_damp * d_shear_condition_dx * d_shear_condition_dx.transpose();

    remove_small_values(shear_force);
    remove_small_values(shear_damping_force);
    remove_small_values(d_shear_force_dx);
    remove_small_values(d_shear_damping_force_dx);
    remove_small_values(d_shear_damping_force_dv);
}
void compute_bend_force(Eigen::Vector12d &bend_force,
                        Eigen::Matrix1212d &d_bend_force_dx, 
                        Eigen::Vector12d &bend_damping_force,
                        Eigen::Matrix1212d &d_bend_damping_force_dx,
                        Eigen::Matrix1212d &d_bend_damping_force_dv,
                        Eigen::Ref<const Eigen::Vector9d> triangle_dot1, 
                        Eigen::Ref<const Eigen::Vector9d> triangle_dot2, 
                        Eigen::Ref<const Eigen::Vector9d> ordered_triangle1, 
                        Eigen::Ref<const Eigen::Vector9d> ordered_triangle2, 
                        double k, double k_damp) {
    double bend_condition;
    Eigen::Vector12d d_bend_condition_dx;
    Eigen::Matrix1212d d2_bend_condition_dx2;
    compute_bend_condition(bend_condition, ordered_triangle1, ordered_triangle2);
    compute_d_bend_condition_dx(d_bend_condition_dx, ordered_triangle1, ordered_triangle2);
    compute_d2_bend_condition_dx2(d2_bend_condition_dx2,
                                    ordered_triangle1, 
                                    ordered_triangle2);
    bend_force = k * bend_condition * d_bend_condition_dx;
    d_bend_force_dx = k * (
                d_bend_condition_dx * d_bend_condition_dx.transpose()
                + d2_bend_condition_dx2 * bend_condition
            );
    Eigen::Vector12d triangle_dot;
    for (int i = 0; i < 12; i++) {
        if (i < 9) {
            triangle_dot(i) = triangle_dot1(i);
        } else {
            triangle_dot(i) = triangle_dot2(i - 3);
        }
    }

   double cdot = d_bend_condition_dx.transpose() * triangle_dot;
    bend_damping_force = k_damp * d_bend_condition_dx.transpose() * cdot;
    d_bend_damping_force_dx = k_damp * d2_bend_condition_dx2 * cdot;
    d_bend_damping_force_dv = k_damp * d_bend_condition_dx * d_bend_condition_dx.transpose();

    remove_small_values(bend_force);
    remove_small_values(bend_damping_force);
    remove_small_values(d_bend_force_dx);
    remove_small_values(d_bend_damping_force_dx);
    remove_small_values(d_bend_damping_force_dv);
}