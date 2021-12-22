#include <compute_bend_condition.h>
#include<iostream>

void triangle_norm(Eigen::Vector3d &norm, Eigen::Ref<const Eigen::Vector9d> triangle_vertices) {
    Eigen::Vector3d xi = triangle_vertices.segment<3>(0);
    Eigen::Vector3d xj = triangle_vertices.segment<3>(3);
    Eigen::Vector3d xk = triangle_vertices.segment<3>(6);         

    norm = (xj - xi).cross(xk-xi).normalized();
}

void compute_bend_condition(double &bend_condition, 
                            Eigen::Ref<const Eigen::Vector9d> triangle_vertices_1, 
                            Eigen::Ref<const Eigen::Vector9d> triangle_vertices_2) {
    //     t02 
    //    /   \
    // t00 - - t01
    // t11 - - t10
    //    \   /
    //     t12
    int common_vertex1 = 0;
    int common_vertex2 = 1;
    // let a
    Eigen::Vector3d n1;
    Eigen::Vector3d n2;
    triangle_norm(n1, triangle_vertices_1);
    triangle_norm(n2, triangle_vertices_2);

    Eigen::Vector3d e = triangle_vertices_1.segment<3>(common_vertex2*3) - triangle_vertices_1.segment<3>(common_vertex1*3);
    e = e.normalized();

    double sin_theta = n1.cross(n2).dot(e);
    double cos_theta = n1.dot(n2);

    bend_condition = atan2(sin_theta, cos_theta);
};

void d_triangle_norm_dx(Eigen::Matrix<double, 3, 12> &result, 
                        Eigen::Ref<const Eigen::Vector9d> triangle_vertices) {

    Eigen::Vector3d xi = triangle_vertices.segment<3>(0);
    Eigen::Vector3d xj = triangle_vertices.segment<3>(3);
    Eigen::Vector3d xk = triangle_vertices.segment<3>(6);         

    Eigen::Matrix39d d_xji_dx;
    d_xji_dx << -1,  0,  0, 1, 0, 0, 0, 0, 0,
                 0, -1,  0, 0, 1, 0, 0, 0, 0,
                 0,  0, -1, 0, 0, 1, 0, 0, 0;

    Eigen::Matrix39d d_xki_dx;
    d_xki_dx << -1,  0,  0, 0, 0, 0, 1, 0, 0,
                 0, -1,  0, 0, 0, 0, 0, 1, 0,
                 0,  0, -1, 0, 0, 0, 0, 0, 1;

    Eigen::Vector3d xji = xj - xi;
    Eigen::Vector3d xki = xk - xi;
    // One of the assumption is to treat this as constant. 
    // So no need to take the derivative for this
    double norm = (xji).cross(xki).norm();

    result.setZero();
    for (int i = 0; i < 9; i++) {
        Eigen::Vector3d tmp1 = d_xji_dx.col(i);
        Eigen::Vector3d temp2 = d_xki_dx.col(i);
        result.col(i) = (tmp1.cross(xki) + xji.cross(temp2))/norm;
    }
}

void d_triangle2_norm_dx(Eigen::Matrix<double, 3, 12> &result, 
                         Eigen::Ref<const Eigen::Vector9d> triangle_vertices) {
    // x here is the common edge in the lower triangle, which are point 0 and point 1
    //     t02 
    //    /   \
    // t00 - - t01
    // t11 - - t10
    //    \   /
    //     t12

    Eigen::Vector3d xi = triangle_vertices.segment<3>(0);
    Eigen::Vector3d xj = triangle_vertices.segment<3>(3);
    Eigen::Vector3d xk = triangle_vertices.segment<3>(6);         
    
    Eigen::Matrix39d d_xji_dx;
    d_xji_dx << 1, 0, 0, -1,  0,  0, 0, 0, 0,
                0, 1, 0,  0, -1,  0, 0, 0, 0,
                0, 0, 1,  0,  0, -1, 0, 0, 0;

    Eigen::Matrix39d d_xki_dx;
    d_xki_dx << 0, 0, 0, -1,  0,  0, 1, 0, 0,
                0, 0, 0,  0, -1,  0, 0, 1, 0,
                0, 0, 0,  0,  0, -1, 0, 0, 1;

    Eigen::Vector3d xji = xj - xi;
    Eigen::Vector3d xki = xk - xi;
    // One of the assumption is to treat this as constant. 
    // So no need to take the derivative for this
    double norm = (xji).cross(xki).norm();
    
    result.setZero();
    for (int i = 0; i < 9; i++) {
        Eigen::Vector3d tmp1 = d_xji_dx.col(i);
        Eigen::Vector3d temp2 = d_xki_dx.col(i);
        if (i >= 6) {
            result.col(i + 3) = (tmp1.cross(xki) + xji.cross(temp2))/norm;
        } else {
            result.col(i) = (tmp1.cross(xki) + xji.cross(temp2))/norm;
        }
    }
}

void compute_d_bend_condition_dx(Eigen::Vector12d &d_bend_condition_dx, 
                                 Eigen::Ref<const Eigen::Vector9d> triangle_vertices_1, 
                                 Eigen::Ref<const Eigen::Vector9d> triangle_vertices_2) {
    int common_vertex1 = 0;
    int common_vertex2 = 1;
    // let a
    Eigen::Vector3d n1;
    Eigen::Vector3d n2;
    triangle_norm(n1, triangle_vertices_1);
    triangle_norm(n2, triangle_vertices_2);

    Eigen::Matrix<double, 3, 12> d_n1_dx;
    Eigen::Matrix<double, 3, 12> d_n2_dx;
    d_triangle_norm_dx(d_n1_dx, triangle_vertices_1);
    d_triangle2_norm_dx(d_n2_dx, triangle_vertices_2);
    
    Eigen::Vector3d e = triangle_vertices_1.segment<3>(common_vertex2*3) - triangle_vertices_1.segment<3>(common_vertex1*3);
    double e_norm = e.norm();
    e = e.normalized();
    Eigen::Matrix<double, 3, 12> d_e_dx;
    d_e_dx << -1,  0,  0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
               0, -1,  0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
               0,  0, -1, 0, 0, 1, 0, 0, 0, 0, 0, 0; 
    d_e_dx = d_e_dx / e_norm;

    double sin_theta = n1.cross(n2).dot(e);
    double cos_theta = n1.dot(n2);

    Eigen::Vector12d d_sin_theta_dx;
    Eigen::Vector12d d_cos_theta_dx;
    for(int i = 0; i < 12; i++) {
        d_sin_theta_dx(i) = (d_n1_dx.col(i).cross(n2) + n1.cross(d_n2_dx.col(i))).dot(e)
                             + n1.cross(n2).dot(d_e_dx.col(i));
        d_cos_theta_dx(i) = d_n1_dx.col(i).dot(n2) + n1.dot(d_n2_dx.col(i));
    }

    double sin2_cos2 = sin_theta * sin_theta + cos_theta * cos_theta;
    double d_atan2_dsin = cos_theta / (sin2_cos2);
    double d_atan2_dcos = - sin_theta / (sin2_cos2);

    d_bend_condition_dx = d_atan2_dcos * d_cos_theta_dx + d_atan2_dsin * d_sin_theta_dx;
}

Eigen::Matrix<double, 36, 12> d2_triangle_norm_dx2_cache;
bool d2_triangle_norm_dx2_cache_cached = false;
void d2_triangle_norm_dx2(Eigen::Matrix<double, 36, 12> &result, 
                          Eigen::Ref<const Eigen::Vector9d> triangle_vertices) {
    Eigen::Vector3d xi = triangle_vertices.segment<3>(0);
    Eigen::Vector3d xj = triangle_vertices.segment<3>(3);
    Eigen::Vector3d xk = triangle_vertices.segment<3>(6);         

    Eigen::Vector3d xji = xj - xi;
    Eigen::Vector3d xki = xk - xi;
    // One of the assumption is to treat this as constant. 
    // So no need to take the derivative for this
    double norm = (xji).cross(xki).norm();
    if (d2_triangle_norm_dx2_cache_cached) {
        result = d2_triangle_norm_dx2_cache / norm;
        return;
    }

    d2_triangle_norm_dx2_cache_cached = true;
    d2_triangle_norm_dx2_cache.setZero();
    Eigen::Matrix39d d_xji_dx;
    d_xji_dx << -1,  0,  0, 1, 0, 0, 0, 0, 0,
                 0, -1,  0, 0, 1, 0, 0, 0, 0,
                 0,  0, -1, 0, 0, 1, 0, 0, 0;

    Eigen::Matrix39d d_xki_dx;
    d_xki_dx << -1,  0,  0, 0, 0, 0, 1, 0, 0,
                 0, -1,  0, 0, 0, 0, 0, 1, 0,
                 0,  0, -1, 0, 0, 0, 0, 0, 1;

    for (int i = 0; i < 9; i++) {
        Eigen::Vector3d temp1 = d_xji_dx.col(i);
        Eigen::Vector3d temp2 = d_xki_dx.col(i);
        for (int j = 0; j < 9; j++) {
            Eigen::Vector3d temp3 = d_xji_dx.col(j);
            Eigen::Vector3d temp4 = d_xki_dx.col(j);
            auto r = (temp1.cross(temp4) + temp3.cross(temp2));
            d2_triangle_norm_dx2_cache(i, j) = r(0);
            d2_triangle_norm_dx2_cache(12 + i, j) = r(1);
            d2_triangle_norm_dx2_cache(24 + i, j) = r(2);
        }
    }
    result = d2_triangle_norm_dx2_cache / norm;
}


Eigen::Matrix<double, 36, 12> d2_triangle2_norm_dx2_cache;
bool d2_triangle2_norm_dx2_cache_cached = false;
void d2_triangle2_norm_dx2(Eigen::Matrix<double, 36, 12> &result, 
                           Eigen::Ref<const Eigen::Vector9d> triangle_vertices) {
    Eigen::Vector3d xi = triangle_vertices.segment<3>(0);
    Eigen::Vector3d xj = triangle_vertices.segment<3>(3);
    Eigen::Vector3d xk = triangle_vertices.segment<3>(6);         

    Eigen::Vector3d xji = xj - xi;
    Eigen::Vector3d xki = xk - xi;
    // One of the assumption is to treat this as constant. 
    // So no need to take the derivative for this
    double norm = (xji).cross(xki).norm();
    if (d2_triangle2_norm_dx2_cache_cached) {
        result = d2_triangle2_norm_dx2_cache / norm;
        return;
    }

    d2_triangle2_norm_dx2_cache_cached = true;
    d2_triangle2_norm_dx2_cache.setZero();
    Eigen::Matrix39d d_xji_dx;
    d_xji_dx << 1, 0, 0, -1,  0,  0, 0, 0, 0,
                0, 1, 0,  0, -1,  0, 0, 0, 0,
                0, 0, 1,  0,  0, -1, 0, 0, 0;

    Eigen::Matrix39d d_xki_dx;
    d_xki_dx << 0, 0, 0, -1,  0,  0, 1, 0, 0,
                0, 0, 0,  0, -1,  0, 0, 1, 0,
                0, 0, 0,  0,  0, -1, 0, 0, 1;
    for (int i = 0; i < 9; i++) {
        Eigen::Vector3d temp1 = d_xji_dx.col(i);
        Eigen::Vector3d temp2 = d_xki_dx.col(i);
        for (int j = 0; j < 9; j++) {
            Eigen::Vector3d temp3 = d_xji_dx.col(j);
            Eigen::Vector3d temp4 = d_xki_dx.col(j);
            auto r = (temp1.cross(temp4) + temp3.cross(temp2));

            int updated_i = (i >= 6) ? i + 3 : i;
            int updated_j = (j >= 6) ? j + 3 : j;
            d2_triangle2_norm_dx2_cache(updated_i, updated_j) = r(0);
            d2_triangle2_norm_dx2_cache(12 + updated_i, updated_j) = r(1);
            d2_triangle2_norm_dx2_cache(24 + updated_i, updated_j) = r(2);
        }
    }
    result = d2_triangle2_norm_dx2_cache / norm;
}

void compute_d2_bend_condition_dx2(Eigen::Matrix1212d &d2_bend_condition_dx2,
                                   Eigen::Ref<const Eigen::Vector9d> triangle_vertices_1, 
                                   Eigen::Ref<const Eigen::Vector9d> triangle_vertices_2) {
    int common_vertex1 = 0;
    int common_vertex2 = 1;
    // let a
    Eigen::Vector3d n1;
    Eigen::Vector3d n2;
    triangle_norm(n1, triangle_vertices_1);
    triangle_norm(n2, triangle_vertices_2);

    Eigen::Matrix<double, 3, 12> d_n1_dx;
    Eigen::Matrix<double, 3, 12> d_n2_dx;
    d_triangle_norm_dx(d_n1_dx, triangle_vertices_1);
    d_triangle2_norm_dx(d_n2_dx, triangle_vertices_2);

    Eigen::Matrix<double, 36, 12> d2_n1_dx2;
    Eigen::Matrix<double, 36, 12> d2_n2_dx2;
    d2_triangle_norm_dx2(d2_n1_dx2, triangle_vertices_1);
    d2_triangle2_norm_dx2(d2_n2_dx2, triangle_vertices_2);
 
    Eigen::Vector3d e = triangle_vertices_1.segment<3>(common_vertex2*3) - triangle_vertices_1.segment<3>(common_vertex1*3);
    double e_norm = e.norm();
    e = e.normalized();

    Eigen::Matrix<double, 3, 12> d_e_dx;
    d_e_dx << -1,  0,  0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
               0, -1,  0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
               0,  0, -1, 0, 0, 1, 0, 0, 0, 0, 0, 0; 
    d_e_dx = d_e_dx / e_norm;

    double sin_theta = n1.cross(n2).dot(e);
    double cos_theta = n1.dot(n2);

    Eigen::Vector12d d_sin_theta_dx;
    Eigen::Vector12d d_cos_theta_dx;

    for(int i = 0; i < 12; i++) {
        d_sin_theta_dx(i) = (d_n1_dx.col(i).cross(n2) + n1.cross(d_n2_dx.col(i))).dot(e)
                             + n1.cross(n2).dot(d_e_dx.col(i));
        d_cos_theta_dx(i) = d_n1_dx.col(i).dot(n2) + n1.dot(d_n2_dx.col(i));
    }

    Eigen::Matrix1212d d2_sin_theta_dx2;
    Eigen::Matrix1212d d2_cos_theta_dx2;
    d2_cos_theta_dx2.setZero();
    d2_sin_theta_dx2.setZero();

    Eigen::Vector3d local_d_n1_dx_i; 
    Eigen::Vector3d local_d_n2_dx_i; 
    Eigen::Vector3d local_d_e_dx_i; 
    Eigen::Vector3d local_d_n1_dx_j;
    Eigen::Vector3d local_d_n2_dx_j; 
    Eigen::Vector3d local_d_e_dx_j; 
    Eigen::Vector3d local_d2_n1_dx2; 
    Eigen::Vector3d local_d2_n2_dx2; 
    for (int i = 0; i < 12; i++) {
        for (int j = i; j < 12; j++) {
            local_d_n1_dx_i = d_n1_dx.col(i);
            local_d_n2_dx_i = d_n2_dx.col(i);
            local_d_e_dx_i = d_e_dx.col(i);

            local_d_n1_dx_j = d_n1_dx.col(j);
            local_d_n2_dx_j = d_n2_dx.col(j);
            local_d_e_dx_j = d_e_dx.col(j);

            local_d2_n1_dx2 << d2_n1_dx2(i, j), d2_n1_dx2(12 + i, j), d2_n1_dx2(24 + i, j);
            local_d2_n2_dx2 << d2_n2_dx2(i, j), d2_n2_dx2(12 + i, j), d2_n2_dx2(24 + i, j);
             
            double local_cos = local_d2_n1_dx2.dot(n2) + n1.dot(local_d2_n2_dx2)
                                + local_d_n1_dx_j.dot(local_d_n2_dx_i) 
                                + local_d_n1_dx_i.dot(local_d_n2_dx_j);

            double local_sin = 
                (
                  local_d2_n1_dx2.cross(n2) + 
                  n1.cross(local_d2_n2_dx2) + 
                  local_d_n1_dx_i.cross(local_d_n2_dx_j) +
                  local_d_n1_dx_j.cross(local_d_n2_dx_i)
                ).dot(e) + 
                (
                  local_d_n1_dx_i.cross(n2) + 
                  n1.cross(local_d_n2_dx_i)
                ).dot(local_d_e_dx_j) +
                (
                  local_d_n1_dx_j.cross(n2) + 
                  n1.cross(local_d_n2_dx_j)
                ).dot(local_d_e_dx_i);

            d2_cos_theta_dx2(i, j) = local_cos;
            d2_sin_theta_dx2(i, j) = local_sin; 
            if (i != j) {
                d2_cos_theta_dx2(j, i) = local_cos;
                d2_sin_theta_dx2(j, i) = local_sin;
            }
        }
    }

    double sin2_cos2 = sin_theta * sin_theta + cos_theta * cos_theta;

    double d_sin_dx_i, d_sin_dx_j, d_cos_dx_i, d_cos_dx_j, d2_sin_dx2_ij, d2_cos_dx2_ij;

    for (int i = 0; i < 12; i++) {
        for (int j = i; j < 12; j++) { 
            d_sin_dx_i = d_sin_theta_dx(i);
            d_sin_dx_j = d_sin_theta_dx(j);
            d_cos_dx_i = d_cos_theta_dx(i);
            d_cos_dx_j = d_cos_theta_dx(j);
            d2_sin_dx2_ij = d2_sin_theta_dx2(i, j);
            d2_cos_dx2_ij = d2_cos_theta_dx2(i, j);
            
            double out = 
                - 2.0 / (sin2_cos2 * sin2_cos2) * 
                ( sin_theta * d_sin_dx_j + cos_theta * d_cos_dx_j) *
                (cos_theta * d_sin_dx_i - sin_theta * d_cos_dx_i)
               +
                1.0 / (sin2_cos2) * 
                (
                    d_cos_dx_j * d_sin_dx_i + cos_theta * d2_sin_dx2_ij -
                    d_sin_dx_j * d_cos_dx_i - sin_theta * d2_cos_dx2_ij
                );

            d2_bend_condition_dx2(i, j) = out;
            if (i != j) {
                d2_bend_condition_dx2(j, i) = out;
            }
        }
    }
}