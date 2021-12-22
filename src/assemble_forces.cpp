#include <assemble_forces.h>
#include <compute_forces.h>
#include <compute_stretch_condition.h>
#include <compute_shear_condition.h>
#include <compute_bend_condition.h>
#include <constants.h>
#include <iostream>

extern double K_BEND;
extern double K_STRETCH;
extern double K_SHEAR;
extern double K_DAMP;
extern double bu;
extern double bv;
void get_triangle(Eigen::Vector9d &triangle, 
                  Eigen::Ref<const Eigen::VectorXd> q,
                  Eigen::Ref<Eigen::RowVectorXi> element) {
    triangle.segment<3>(0) = q.segment<3>(element(0) * 3);
    triangle.segment<3>(3) = q.segment<3>(element(1) * 3);
    triangle.segment<3>(6) = q.segment<3>(element(2) * 3); 
}

void order_triangles(Eigen::Vector3i &ordered_element1, 
                     Eigen::Vector3i &ordered_element2,
                     Eigen::Ref<const Eigen::Vector3i> element1,
                     Eigen::Ref<const Eigen::Vector3i> element2) {
    //     e02 
    //    /   \
    // e00 - - e01
    // e11 - - e10
    //    \   /
    //     e12
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            if (element1(i) != element2(j)) {
                continue;
            }
            if (element1((i + 1)%3) == element2((j + 2) % 3)) {
                ordered_element1(0) = element1(i);
                ordered_element1(1) = element1((i +1)%3);
                ordered_element1(2) = element1((i+2)%3);
                ordered_element2(0) = element2((j+2)%3);
                ordered_element2(1) = element2(j);
                ordered_element2(2) = element2((j+1)%3);
            } else {
                ordered_element1(0) = element1((i+2)%3);
                ordered_element1(1) = element1(i);
                ordered_element1(2) = element1((i+1)%3);
                ordered_element2(0) = element2(j);
                ordered_element2(1) = element2((j+1)%3);
                ordered_element2(2) = element2((j+2)%3);
            }
            return;
        }
    }
}
void assemble_forces(Eigen::VectorXd &f, 
                     Eigen::SparseMatrixd &df_dx, 
                     Eigen::SparseMatrixd &df_dv, 
                     Eigen::Ref<const Eigen::VectorXd> q,
                     Eigen::Ref<const Eigen::VectorXd> qdot,
                     Eigen::Ref<const Eigen::MatrixXi> F,
                     Eigen::Ref<const Eigen::MatrixXd> uv_coords,
                     const std::vector<std::vector<int>> &adjacency_list,
                     Eigen::Ref<const Eigen::VectorXd> a0) { 
    f.resize(q.size());  
    f.setZero();  
    df_dx.resize(q.size(), q.size());  
    df_dx.setZero();  
    df_dv.resize(q.size(), q.size());  
    df_dv.setZero();  
    std::vector<Eigen::Triplet<double>> update_list;
    std::vector<Eigen::Triplet<double>> update_list_dv;
    for (int triangle_index = 0; triangle_index < F.rows(); triangle_index++) {

        Eigen::RowVectorXi element = F.row(triangle_index);

        Eigen::Vector9d triangle;
        get_triangle(triangle, q, element);
        Eigen::Vector9d triangle_dot;
        get_triangle(triangle_dot, qdot, element);

        Eigen::Vector6d ref_plane_coordinate = uv_coords.row(triangle_index);
        double area = a0(triangle_index);
        // common
        Eigen::Matrix32d stretch;
        compute_stretch(stretch, ref_plane_coordinate, triangle);
        Eigen::Matrix39d d_stretch_u_dx;
        Eigen::Matrix39d d_stretch_v_dx;
        d_stretch_dx(d_stretch_u_dx, d_stretch_v_dx, ref_plane_coordinate);

        //stretch
        Eigen::Vector9d stretch_force;
        Eigen::Matrix99d d_stretch_force_dx;
        Eigen::Vector9d stretch_damping_force;
        Eigen::Matrix99d d_stretch_damping_force_dx;
        Eigen::Matrix99d d_stretch_damping_force_dv;

        compute_stretch_force(stretch_force, d_stretch_force_dx,
                              stretch_damping_force, d_stretch_damping_force_dx, d_stretch_damping_force_dv,
                              triangle_dot, stretch, d_stretch_u_dx, d_stretch_v_dx,
                              K_STRETCH, (K_DAMP * K_STRETCH), area, bu, bv);

        // shear
        Eigen::Vector9d shear_force;
        Eigen::Matrix99d d_shear_force_dx;
        Eigen::Vector9d shear_damping_force ;
        Eigen::Matrix99d d_shear_damping_force_dx;
        Eigen::Matrix99d d_shear_damping_force_dv;
        compute_shear_force(shear_force, d_shear_force_dx, shear_damping_force, d_shear_damping_force_dx, d_shear_damping_force_dv,
                            triangle_dot, stretch, d_stretch_u_dx, d_stretch_v_dx,
                            K_SHEAR, (K_DAMP * K_SHEAR), area);

        for (int i = 0; i < 9; i++) {
            int element_i = element(i/3) * 3 + i%3;
            f(element_i) -= (stretch_force(i) +
                            stretch_damping_force(i) +
                            shear_force(i) +
                            shear_damping_force(i));
        }
            // if (i % 3 == 0 && f(element_i) != 0) {
                // std::cout << "WHY???";
            // }
        Eigen::Matrix99d sum_dforce = d_stretch_force_dx
                                  + d_stretch_damping_force_dx
                                  + d_shear_force_dx
                                  + d_shear_damping_force_dx;
         Eigen::Matrix99d sum_dv = d_stretch_damping_force_dv
                                  + d_shear_damping_force_dv;

        for (int i = 0; i < 3; i++) {
            int index_i = element(i) * 3;
            for (int j = 0; j < 3; j++) {
                int index_j = element(j) * 3;
                for (int k = 0; k < 3; k++) {
                    update_list.push_back(Eigen::Triplet<double>(index_i, index_j + k, -sum_dforce(i * 3, j * 3 + k)));
                    update_list.push_back(Eigen::Triplet<double>(index_i + 1, index_j + k, -sum_dforce(i * 3 + 1, j * 3 + k)));
                    update_list.push_back(Eigen::Triplet<double>(index_i + 2, index_j + k, -sum_dforce(i * 3 + 2, j * 3 + k)));
 
                    update_list_dv.push_back(Eigen::Triplet<double>(index_i, index_j + k, -sum_dv(i * 3, j * 3 + k)));
                    update_list_dv.push_back(Eigen::Triplet<double>(index_i + 1, index_j + k, -sum_dv(i * 3 + 1, j * 3 + k)));
                    update_list_dv.push_back(Eigen::Triplet<double>(index_i + 2, index_j + k, -sum_dv(i * 3 + 2, j * 3 + k)));
                }
            }
        }

        // bend
        std::vector<int> adj_list = adjacency_list.at(triangle_index);
        for (int i = 0; i < adj_list.size(); i++) {
            if (adj_list.at(i) < triangle_index) {
                // avoid double count;
                exit(1);
                // continue;
            }
            Eigen::RowVectorXi element2 = F.row(adj_list.at(i));
            Eigen::Vector3i ordered_element1;
            Eigen::Vector3i ordered_element2;
            order_triangles(ordered_element1, ordered_element2, element, element2);

            Eigen::Vector9d ordered_triangle1;
            Eigen::Vector9d ordered_triangle2;
            get_triangle(ordered_triangle1, q, ordered_element1);
            get_triangle(ordered_triangle2, q, ordered_element2);

            Eigen::Vector9d triangle_dot1;
            Eigen::Vector9d triangle_dot2;
            get_triangle(triangle_dot1, qdot, ordered_element1);
            get_triangle(triangle_dot2, qdot, ordered_element2);

            Eigen::Vector12d bend_force;
            Eigen::Matrix1212d d_bend_force_dx;
            Eigen::Vector12d bend_damping_force;
            Eigen::Matrix1212d d_bend_damping_force_dx;
            Eigen::Matrix1212d d_bend_damping_force_dv;
            compute_bend_force(bend_force, d_bend_force_dx, bend_damping_force, d_bend_damping_force_dx, d_bend_damping_force_dv,
                        triangle_dot1, triangle_dot2, ordered_triangle1, ordered_triangle2,
                        K_BEND, (K_BEND * K_DAMP));

            for (int x = 0; x < 12; x++) {
                int element_x = x < 9 ? ordered_element1(x/3) : ordered_element2(2);
                element_x = element_x * 3 + x%3;

                f(element_x) -= (bend_force(x) + bend_damping_force(x));
                for (int y = 0; y < 12; y++) {
                    int element_y = y < 9 ? ordered_element1(y/3) : ordered_element2(2);
                    element_y = element_y * 3 + y%3;

                    double dx_value = d_bend_force_dx(x, y) + d_bend_damping_force_dx(x, y);
                    double dv_value = d_bend_damping_force_dv(x, y);
                    update_list.push_back(Eigen::Triplet<double>(element_x, element_y, -dx_value));
                    update_list_dv.push_back(Eigen::Triplet<double>(element_x, element_y, -dv_value));
                }
            }
        }
    }

    df_dx.setFromTriplets(update_list.begin(), update_list.end());   
    df_dv.setFromTriplets(update_list_dv.begin(), update_list_dv.end());   
};
