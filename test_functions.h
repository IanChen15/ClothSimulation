#ifndef TEST_FUNCTIONS_H
#define TEST_FUNCTIONS_H

//Assignment 4 code 
#include <igl/readMESH.h>
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/readOFF.h>
#include <read_tetgen.h>
#include <igl/boundary_facets.h>
#include <igl/volume.h>

#include <unit_test.h>

inline void test_functions() {
    unit_test();
}

// void test_stretch_forces() {
//     std::cout << "TEST stretch_forces \n";
//     Eigen::Vector6d ref_plane_coordinate;
//     ref_plane_coordinate << 0, 0, 4, 0, 2, 2;

//     Eigen::Vector9d triangle_vertices;
//     triangle_vertices << 0, 0, 3, 4, 0, 3, 2, 2, 3;

//     double area = 4;

//     Eigen::Matrix32d output;
//     stretch_forces(output, ref_plane_coordinate, triangle_vertices);
//     Eigen::Matrix32d expected;
//     expected << 1, 0, 0, 1, 0, 0;
 
//     assertTrue(expected == output);
// }

// void test_compute_stretch_condition() {
//     std::cout << "TEST compute_stretch_condition \n";

//     double area = 4;
//     double bu = 0.2;
//     double bv = 0.3;
//     Eigen::Matrix32d stretch_force_matrix;
//     stretch_force_matrix << 1, 0, 0, 1, 0, 0;

//     Eigen::Vector2d output;
//     compute_stretch_condition(output, stretch_force_matrix, area, bu, bv);
//     // compute_stretch_condition();

//     Eigen::Vector2d expected;
//     expected << 3.2, 2.8;  
//     assertTrue(expected == output);
// }

// void test_compute_shear_condition() {
//     std::cout << "TEST compute_shear_condition \n";

//     double area = 2;
//     Eigen::Matrix32d stretch_force_matrix;
//     stretch_force_matrix << 1, 2, 3, 4, 5, 6;

//     double output;
//     compute_shear_condition(output, stretch_force_matrix, area);

//     double expected = 88;
//     assertTrue(expected == output);
// }
#endif

