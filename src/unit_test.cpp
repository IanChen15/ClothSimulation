#include <unit_test.h>
#include<iostream>
#include <compute_stretch_condition.h>
#include <compute_shear_condition.h>
#include <compute_bend_condition.h>
#include <assemble_forces.h>
#include <compute_forces.h>

void test_stretch_forces();
void test_compute_stretch_condition();
void test_compute_shear_condition();
void test_compute_bend_condition();
void test_d_stretch_dx();
void test_d_stretch_conditiondx();
void test_d_shear_condition_dx();
void test_d_triangle_norm_dx();
void test_d_triangle2_norm_dx();
void test_compute_d_bend_condition_dx();
void test_compute_d2_stretch_condition_dx2();
void test_d2_shear_condition_dx2();
void test_d2_triangle_norm_dx2();
void test_compute_d2_bend_condition_dx2();
void test_order_triangles();
void test_shear_force();
void test_bend_force();
void test_stret_force();
int test_count = 1;
void assertTrue(bool b) {
    if (b)
        std::cout << "    Test " << test_count++ << ": True\n";
    else {
        std::cout << "    Test " << test_count++ << ": Failed \n";
        exit(1);
    }
}

void unit_test() {
    test_stretch_forces();
    test_compute_stretch_condition();
    test_compute_shear_condition();
    test_compute_bend_condition();

    test_d_stretch_dx();
    test_d_stretch_conditiondx();
    test_d_shear_condition_dx();

    test_d_triangle_norm_dx();
    test_d_triangle2_norm_dx();

    test_compute_d_bend_condition_dx();
    test_compute_d2_stretch_condition_dx2();
    test_d2_shear_condition_dx2();
    test_d2_triangle_norm_dx2();
    // test_d2_triangle2_norm_dx2();
    test_compute_d2_bend_condition_dx2();
    
    test_order_triangles();
    test_shear_force();
    test_bend_force();
    test_stret_force();
}

void test_stretch_forces() {
    std::cout << "TEST stretch_forces \n";
    Eigen::Vector6d ref_plane_coordinate;
    ref_plane_coordinate << 0, 0, 4, 0, 2, 2;

    Eigen::Vector9d triangle_vertices;
    triangle_vertices << 0, 0, 3, 4, 0, 3, 2, 2, 3;

    double area = 4;

    Eigen::Matrix32d output;
    compute_stretch(output, ref_plane_coordinate, triangle_vertices);
    Eigen::Matrix32d expected;
    expected << 1, 0, 0, 1, 0, 0;
 
    assertTrue(expected == output);
}

void test_compute_stretch_condition() {
    std::cout << "TEST compute_stretch_condition \n";

    double area = 4;
    double bu = 0.2;
    double bv = 0.3;
    Eigen::Matrix32d stretch_force_matrix;
    stretch_force_matrix << 1, 0, 0, 1, 0, 0;

    Eigen::Vector2d output;
    compute_stretch_condition(output, stretch_force_matrix, area, bu, bv);
    // compute_stretch_condition();

    Eigen::Vector2d expected;
    expected << 3.2, 2.8;  
    assertTrue(expected == output);
}


void test_compute_shear_condition() {
    std::cout << "TEST compute_shear_condition \n";

    double area = 2;
    Eigen::Matrix32d stretch_force_matrix;
    stretch_force_matrix << 1, 2, 3, 4, 5, 6;

    double output;
    compute_shear_condition(output, stretch_force_matrix, area);

    double expected = 88;
    assertTrue(expected == output);
}

void test_compute_bend_condition() {
    std::cout << "TEST compute_bend_condition \n";

    Eigen::Vector9d triangle_vertices1;
    Eigen::Vector9d triangle_vertices2;
    triangle_vertices1 << 0, 0, 3, 4, 0, 3, 2, 2, 3;
    triangle_vertices2 << 4, 0, 3, 0, 0, 3, 2, 0, 1;

    double output;
    compute_bend_condition(output, triangle_vertices1, triangle_vertices2);

    double expected = atan2(1, 0);
    assertTrue(expected == output);
}


void test_d_stretch_dx() {
    std::cout << "TEST d_stretch_dx \n";
    Eigen::Vector6d ref_plane_coordinate;
    ref_plane_coordinate << 0, 0, 4, 0, 2, 2;

    Eigen::Matrix39d output1;
    Eigen::Matrix39d output2;
    d_stretch_dx(output1, output2, ref_plane_coordinate); 

    Eigen::Matrix39d expected1;
    expected1 << -0.25,     0,     0, 0.25,    0,    0, 0, 0, 0,
                     0, -0.25,     0,    0, 0.25,    0, 0, 0, 0,
                     0,     0, -0.25,    0,    0, 0.25, 0, 0, 0;

    Eigen::Matrix39d expected2;
    expected2 << -0.25,     0,     0, -0.25,     0,     0, 0.5,   0,   0,
                     0, -0.25,     0,     0, -0.25,     0,   0, 0.5,   0,
                     0,     0, -0.25,     0,     0, -0.25,   0,   0, 0.5;
    assertTrue(expected1 == output1);
    assertTrue(expected2 == output2);
}

void test_d_stretch_conditiondx() {
    std::cout << "TEST d_stretch_condition_dx \n";
    Eigen::Vector6d ref_plane_coordinate;
    ref_plane_coordinate << 0, 0, 4, 0, 2, 2;
    Eigen::Vector9d triangle_vertices;
    triangle_vertices << 0, 0, 3, 4, 0, 3, 2, 2, 3;
    double area = 4;

    Eigen::Matrix39d d_stretch_u_dx;
    Eigen::Matrix39d d_stretch_v_dx;
    d_stretch_dx(d_stretch_u_dx, d_stretch_v_dx, ref_plane_coordinate); 

    Eigen::Matrix32d stretch_force;
    compute_stretch(stretch_force, ref_plane_coordinate, triangle_vertices);

    Eigen::Matrix<double, 2, 9> output;
    compute_d_stretch_condition_dx(output, 
                                   d_stretch_u_dx, 
                                   d_stretch_v_dx,
                                   stretch_force, 
                                   area);

    Eigen::Matrix<double, 2, 9> expected;
    expected << -1,  0, 0, 1,  0, 0, 0, 0, 0,
                 0, -1, 0, 0, -1, 0, 0, 2, 0;

    assertTrue(expected == output);
}

void test_compute_d2_stretch_condition_dx2() {
    std::cout << "TEST compute_d2_stretch_condition_dx2\n";
    Eigen::Vector6d ref_plane_coordinate;
    ref_plane_coordinate << 0, 0, 4, 0, 2, 2;
    Eigen::Vector9d triangle_vertices;
    triangle_vertices << 0, 0, 3, 4, 0, 3, 2, 2, 3;
    double area = 4;

    Eigen::Matrix39d d_stretch_u_dx;
    Eigen::Matrix39d d_stretch_v_dx;
    d_stretch_dx(d_stretch_u_dx, d_stretch_v_dx, ref_plane_coordinate); 

    Eigen::Matrix32d stretch_force;
    compute_stretch(stretch_force, ref_plane_coordinate, triangle_vertices);
    
    Eigen::Matrix99d d2_stretch_condition_u_dx2; 
    Eigen::Matrix99d d2_stretch_condition_v_dx2;
    
    compute_d2_stretch_condition_dx2(d2_stretch_condition_u_dx2, 
                                     d2_stretch_condition_v_dx2,
                                     d_stretch_u_dx, 
                                     d_stretch_v_dx,
                                     stretch_force, 
                                     area);
}

void test_d_shear_condition_dx() {
    std::cout << "TEST d_stretch_condition_dx \n";
    Eigen::Vector6d ref_plane_coordinate;
    ref_plane_coordinate << 0, 0, 4, 0, 2, 2;
    Eigen::Vector9d triangle_vertices;
    triangle_vertices << 0, 0, 3, 4, 0, 3, 2, 2, 3;
    double area = 4;

    Eigen::Matrix39d d_stretch_u_dx;
    Eigen::Matrix39d d_stretch_v_dx;
    d_stretch_dx(d_stretch_u_dx, d_stretch_v_dx, ref_plane_coordinate); 

    Eigen::Matrix32d stretch_force;
    compute_stretch(stretch_force, ref_plane_coordinate, triangle_vertices);

    Eigen::Vector9d output;
    compute_d_shear_condition_dx(output, 
                                   d_stretch_u_dx, 
                                   d_stretch_v_dx,
                                   stretch_force, 
                                   area);

    Eigen::Vector9d expected;
    expected << -1, -1, 0, -1, 1, 0, 2, 0, 0;
    assertTrue(expected == output);
}

void test_d_triangle_norm_dx() {
    std::cout << "TEST d_triangle_norm_dx \n";

    Eigen::Vector9d triangle_vertices1;
    triangle_vertices1 << 0, 0, 3, 4, 0, 3, 2, 2, 3;
    double norm = 8;

    Eigen::Matrix<double, 3, 12> output;
    d_triangle_norm_dx(output, triangle_vertices1);

    Eigen::Matrix<double, 3, 12> expected;
    expected << 0,  0, 2, 0,  0, -2, 0, 0,  0, 0, 0, 0, 
                0,  0, 2, 0,  0,  2, 0, 0, -4, 0, 0, 0, 
               -2, -2, 0, 2, -2,  0, 0, 4,  0, 0, 0, 0;
    expected /= norm;
    assertTrue(expected == output);
}

void test_d_triangle2_norm_dx() {
    std::cout << "TEST d_triangle2_norm_dx \n";

    Eigen::Vector9d triangle_vertices2;
    triangle_vertices2 << 4, 0, 3, 0, 0, 3, 2, 0, 1;
    double norm = 8;

    Eigen::Matrix<double, 3, 12>  output;
    d_triangle2_norm_dx(output, triangle_vertices2);

   Eigen::Matrix<double, 3, 12>  expected;
    expected << 0, -2,  0,  0, 2,  0, 0, 0, 0, 0,  0, 0,
                2,  0, -2, -2, 0, -2, 0, 0, 0, 0,  0, 4, 
                0,  2,  0,  0, 2,  0, 0, 0, 0, 0, -4, 0;
    expected /= norm;
    assertTrue(expected == output);
}


void test_compute_d_bend_condition_dx() {
    std::cout << "TEST compute_d_bend_condition_dx \n";

    Eigen::Vector9d triangle_vertices1;
    Eigen::Vector9d triangle_vertices2;
    triangle_vertices1 << 1, 1, 3, 4, 0, 3, 2, 2, 3;
    triangle_vertices2 << 4, 0, 3, 1, 1, 3, 2, 0, 1;

    Eigen::Vector12d output;
    compute_d_bend_condition_dx(output, triangle_vertices1, triangle_vertices2);
    // Eigen::Vector12d output2;
    // compute_d_bend_condition_dx(output2, triangle_vertices2, triangle_vertices1);
    // assertTrue(output.segment<6>(0) == -output2.segment<6>(0));
    // assertTrue(output.segment<3>(6) == -output2.segment<3>(9));
    // assertTrue(output.segment<3>(9) == -output2.segment<3>(6));
    std::cout << output << "\n\n";
    // std::cout << output2 << "\n\n";

    Eigen::Vector12d expected;
    // result coppied from matlab
    expected << -0.17248787, -0.22998382,  0.71869946,
                 0.02874797, -0.20123585,  0.21560984, 
                          0,           0, -0.79056941, 
                 0.14373989,  0.43121968, -0.14373989;

    
    // check if they have the same string value. Cannot compare directly because 
    // of rounding precision errors.
    std::stringstream ss;
    ss << expected;
    auto str_expected =  ss.str();

    std::stringstream ss2;
    ss2 << output;
    auto str_output =  ss2.str();

    assertTrue(str_expected == str_output);
}


void test_d2_shear_condition_dx2() {
    std::cout << "TEST test_d2_shear_condition_dx2 \n";
    Eigen::Vector6d ref_plane_coordinate;
    ref_plane_coordinate << 0, 0, 4, 0, 2, 2;
    double area = 4;

    Eigen::Matrix39d d_stretch_u_dx;
    Eigen::Matrix39d d_stretch_v_dx;
    d_stretch_dx(d_stretch_u_dx, d_stretch_v_dx, ref_plane_coordinate); 

    Eigen::Matrix99d output;
    compute_d2_shear_condition_dx2(output, 
                                   d_stretch_u_dx, 
                                   d_stretch_v_dx,
                                   area);

    Eigen::Matrix99d expected;
    expected << 0.5,    0,    0,    0,    0,    0, -0.5,    0,    0,
                  0,  0.5,    0,    0,    0,    0,    0, -0.5,    0,
                  0,    0,  0.5,    0,    0,    0,    0,    0, -0.5,
                  0,    0,    0, -0.5,    0,    0,  0.5,    0,    0,
                  0,    0,    0,    0, -0.5,    0,    0,  0.5,    0,
                  0,    0,    0,    0,    0, -0.5,    0,    0,  0.5,
               -0.5,    0,    0,  0.5,    0,    0,    0,    0,    0,
                  0, -0.5,    0,    0,  0.5,    0,    0,    0,    0,
                  0,    0, -0.5,    0,    0,  0.5,    0,    0,    0;
    assertTrue(expected == output);
}

void test_d2_triangle_norm_dx2() {
    std::cout << "TEST test_d2_triangle_norm_dx2 \n";

    Eigen::Vector9d triangle_vertices1;
    triangle_vertices1 << 0, 0, 3, 4, 0, 3, 2, 2, 3;
    double norm = 8;

    Eigen::Matrix<double, 36, 12> output;
    d2_triangle_norm_dx2(output, triangle_vertices1);
    Eigen::Matrix<double, 36, 12> expected;
    expected << 
     0,      0,      0,      0,      0,      0,      0,      0,      0, 0, 0, 0,
     0,      0,      0,      0,      0,  0.125,      0,      0, -0.125, 0, 0, 0,
     0,      0,      0,      0, -0.125,      0,      0,  0.125,      0, 0, 0, 0,
     0,      0,      0,      0,      0,      0,      0,      0,      0, 0, 0, 0,
     0,      0, -0.125,      0,      0,      0,      0,      0,  0.125, 0, 0, 0,
     0,  0.125,      0,      0,      0,      0,      0, -0.125,      0, 0, 0, 0,
     0,      0,      0,      0,      0,      0,      0,      0,      0, 0, 0, 0,
     0,      0,  0.125,      0,      0, -0.125,      0,      0,      0, 0, 0, 0,
     0, -0.125,      0,      0,  0.125,      0,      0,      0,      0, 0, 0, 0,
     0,      0,      0,      0,      0,      0,      0,      0,      0, 0, 0, 0,
     0,      0,      0,      0,      0,      0,      0,      0,      0, 0, 0, 0,
     0,      0,      0,      0,      0,      0,      0,      0,      0, 0, 0, 0,
     0,      0,      0,      0,      0, -0.125,      0,      0,  0.125, 0, 0, 0,
     0,      0,      0,      0,      0,      0,      0,      0,      0, 0, 0, 0,
     0,      0,      0,  0.125,      0,      0, -0.125,      0,      0, 0, 0, 0,
     0,      0,  0.125,      0,      0,      0,      0,      0, -0.125, 0, 0, 0,
     0,      0,      0,      0,      0,      0,      0,      0,      0, 0, 0, 0,
-0.125,      0,      0,      0,      0,      0,  0.125,      0,      0, 0, 0, 0,
     0,      0, -0.125,      0,      0,  0.125,      0,      0,      0, 0, 0, 0,
     0,      0,      0,      0,      0,      0,      0,      0,      0, 0, 0, 0,
 0.125,      0,      0, -0.125,      0,      0,      0,      0,      0, 0, 0, 0,
     0,      0,      0,      0,      0,      0,      0,      0,      0, 0, 0, 0,
     0,      0,      0,      0,      0,      0,      0,      0,      0, 0, 0, 0,
     0,      0,      0,      0,      0,      0,      0,      0,      0, 0, 0, 0,
     0,      0,      0,      0,  0.125,      0,      0, -0.125,      0, 0, 0, 0,
     0,      0,      0, -0.125,      0,      0,  0.125,      0,      0, 0, 0, 0,
     0,      0,      0,      0,      0,      0,      0,      0,      0, 0, 0, 0,
     0, -0.125,      0,      0,      0,      0,      0,  0.125,      0, 0, 0, 0,
 0.125,      0,      0,      0,      0,      0, -0.125,      0,      0, 0, 0, 0,
     0,      0,      0,      0,      0,      0,      0,      0,      0, 0, 0, 0,
     0,  0.125,      0,      0, -0.125,      0,      0,      0,      0, 0, 0, 0,
-0.125,      0,      0,  0.125,      0,      0,      0,      0,      0, 0, 0, 0,
     0,      0,      0,      0,      0,      0,      0,      0,      0, 0, 0, 0,
     0,      0,      0,      0,      0,      0,      0,      0,      0, 0, 0, 0,
     0,      0,      0,      0,      0,      0,      0,      0,      0, 0, 0, 0,
     0,      0,      0,      0,      0,      0,      0,      0,      0, 0, 0, 0;
    assertTrue(expected == output);
}

void test_compute_d2_bend_condition_dx2() {
    std::cout << "TEST d2_bend_condition_dx2 \n";

    Eigen::Vector9d triangle_vertices1;
    Eigen::Vector9d triangle_vertices2;
    triangle_vertices1 << 1, 1, 3, 4, 0, 3, 2, 2, 3;
    triangle_vertices2 << 4, 0, 3, 1, 1, 3, 2, 0, 1;

    Eigen::Matrix1212d output;
    compute_d2_bend_condition_dx2(output, triangle_vertices1, triangle_vertices2);

    Eigen::Matrix1212d expected;
    expected <<
     -0.13067263058547022033053279109226, -0.078403578351282132198319674655357,   0.56842594304679545843781764125134,   0.06010941006931630135204508390244,  0.010453810446837617626442623287381,  -0.14569998310279929566854406206787,                                   0,                                   0, -0.35216273942784224379078587199364,  0.070563220516153918978487707189821,  0.067949767904444514571877051367976, -0.070563220516153918978487707189821,
    -0.078403578351282132198319674655357,   0.18294168281965830846274590752917,    0.3985515232856841720081250128314, -0.015680715670256426439663934931071, -0.033974883952222257285938525683988,   0.10519146762130352736607889682927,                                   0,                                   0, -0.40965869688544914073622030007424,  0.094084294021538558637983609586428,  -0.14896679886743605117680738184518, -0.094084294021538558637983609586428,
      0.56842594304679545843781764125134,    0.3985515232856841720081250128314, -0.047042147010769279318991804793214, -0.094737657174465909739636273541889,  0.088204025645192398723109633987276, -0.031361431340512852879327869862143, -0.39528470752104741649986169305409, -0.39528470752104741649986169305409,                                   0, -0.078403578351282132198319674655357, -0.091470841409829154231372953764583,  0.078403578351282132198319674655357,
      0.06010941006931630135204508390244, -0.015680715670256426439663934931071, -0.094737657174465909739636273541889,  0.023521073505384639659495902396607,   0.12283227275034200711070082362673, -0.071543265245544945630966703123013,                                   0,                                   0,   0.1545203856673185355408550254666, -0.083630483574700941011540986299047,  -0.10715155708008558067103688869565,  0.011760536752692319829747951198303,
     0.010453810446837617626442623287381, -0.033974883952222257285938525683988,  0.088204025645192398723109633987276,   0.12283227275034200711070082362673, 0.0026134526117094044066106558218452,   0.17738809601977582409869826390774,                                   0,                                   0,  -0.1832683643961219840135722395069,  -0.13328608319717962473714344691411,  0.031361431340512852879327869862143, -0.082323757268846238808235658388124,
     -0.14569998310279929566854406206787,   0.10519146762130352736607889682927, -0.031361431340512852879327869862143, -0.071543265245544945630966703123013,   0.17738809601977582409869826390774, -0.020907620893675235252885246574762,  0.19764235376052370824993084652704, -0.19764235376052370824993084652704,                                   0,  0.019600894587820533049579918663839,  -0.08493720988055564321484631420997,  0.052269052234188088132213116436904,
                                       0,                                    0,  -0.39528470752104741649986169305409,                                    0,                                    0,   0.19764235376052370824993084652704,                                   0,                                   0,  0.19764235376052370824993084652704,                                    0,                                    0,                                    0,
                                       0,                                    0,  -0.39528470752104741649986169305409,                                    0,                                    0,  -0.19764235376052370824993084652704,                                   0,                                   0,  0.59292706128157112474979253958113,                                    0,                                    0,                                    0,
     -0.35216273942784224379078587199364,  -0.40965869688544914073622030007424,                                    0,    0.1545203856673185355408550254666,   -0.1832683643961219840135722395069,                                    0,  0.19764235376052370824993084652704,  0.59292706128157112474979253958113,                                   0,                                    0,                                    0,                                    0,
     0.070563220516153918978487707189821,  0.094084294021538558637983609586428, -0.078403578351282132198319674655357, -0.083630483574700941011540986299047,  -0.13328608319717962473714344691411,  0.019600894587820533049579918663839,                                   0,                                   0,                                   0,  0.013067263058547022033053279109226,  0.039201789175641066099159837327678,  0.058802683763461599148739755991517,
     0.067949767904444514571877051367976,  -0.14896679886743605117680738184518, -0.091470841409829154231372953764583,  -0.10715155708008558067103688869565,  0.031361431340512852879327869862143,  -0.08493720988055564321484631420997,                                   0,                                   0,                                   0,  0.039201789175641066099159837327678,   0.11760536752692319829747951198303,   0.17640805129038479744621926797455,
    -0.070563220516153918978487707189821, -0.094084294021538558637983609586428,  0.078403578351282132198319674655357,  0.011760536752692319829747951198303, -0.082323757268846238808235658388124,  0.052269052234188088132213116436904,                                   0,                                   0,                                   0,  0.058802683763461599148739755991517,   0.17640805129038479744621926797455,  -0.13067263058547022033053279109226;
    bool same = true;
    for(int i = 0; i < 12; i++) {
        for (int j = 0; j< 12; j++) {
            if (abs(expected(i, j) - output(i, j)) > 1e-15) {
                // allow small rounding error
                same = false;
                break;
            }
        }
    }
    assertTrue(same);
}

void test_order_triangles() {
    std::cout << "TEST test_order_triangles\n";
    Eigen::Vector3i element1;
    Eigen::Vector3i element2;
    Eigen::Vector3i ordered_element1;
    Eigen::Vector3i ordered_element2;
    Eigen::Vector3i expected1;
    Eigen::Vector3i expected2;

    element1 << 1, 2, 3;
    element2 << 2, 1, 5;
    expected1 << 1, 2, 3;
    expected2 << 2, 1, 5;
    order_triangles(ordered_element1, 
                    ordered_element2,
                    element1,
                    element2);
    assertTrue(expected1 == ordered_element1 &&
                expected2 == ordered_element2);
    element1 << 1, 2, 3;
    element2 << 5, 2, 1;
    order_triangles(ordered_element1, 
                    ordered_element2,
                    element1,
                    element2);
    assertTrue(expected1 == ordered_element1 &&
                expected2 == ordered_element2);
    element1 << 3, 1, 2;
    element2 << 5, 2, 1;
    order_triangles(ordered_element1, 
                    ordered_element2,
                    element1,
                    element2);
    assertTrue(expected1 == ordered_element1 &&
                expected2 == ordered_element2);
    element1 << 2, 3, 1;
    element2 << 5, 2, 1;
    order_triangles(ordered_element1, 
                    ordered_element2,
                    element1,
                    element2);
    assertTrue(expected1 == ordered_element1 &&
                expected2 == ordered_element2);
}
void test_stret_force() {
    std::cout << "TEST test_stret_force\n";
    Eigen::Vector6d ref_plane_coordinate;
    ref_plane_coordinate << 0, 0, 4, 0, 2, 2;

    Eigen::Vector9d triangle_vertices;
    triangle_vertices << 2, 1, 3, 5, 0, 3, 1, 1, 3;

    double area = 4;

    Eigen::Matrix32d stretch;
    compute_stretch(stretch, ref_plane_coordinate, triangle_vertices);
    Eigen::Matrix39d d_stretch_u_dx;
    Eigen::Matrix39d d_stretch_v_dx;
    d_stretch_dx(d_stretch_u_dx, d_stretch_v_dx, ref_plane_coordinate);

    Eigen::Vector9d stretch_force;
    Eigen::Matrix99d d_stretch_force_dx;
    Eigen::Vector9d stretch_damping_force;
    Eigen::Matrix99d d_stretch_damping_force_dx;
    Eigen::Matrix99d d_stretch_damping_force_dv;
    Eigen::Vector9d triangle_dot;
    triangle_dot <<0, 0, 0, 0, 0, 0, 0, 0, 0;

    compute_stretch_force(stretch_force, 
                          d_stretch_force_dx,
                          stretch_damping_force, 
                          d_stretch_damping_force_dx,
                          d_stretch_damping_force_dv,
                          triangle_dot,
                          stretch,
                          d_stretch_u_dx,
                          d_stretch_v_dx,
                          1, 1, area, 1, 1);
}

void test_shear_force() {
    std::cout << "TEST shear_force\n";
    Eigen::Vector6d ref_plane_coordinate;
    ref_plane_coordinate << 0, 0, 4, 0, 2, 2;

    Eigen::Vector9d triangle_vertices;
    triangle_vertices << 0, 0, 3, 4, 0, 3, 1, 1, 3;

    double area = 4;

    Eigen::Matrix32d stretch;
    compute_stretch(stretch, ref_plane_coordinate, triangle_vertices);
    Eigen::Matrix39d d_stretch_u_dx;
    Eigen::Matrix39d d_stretch_v_dx;
    d_stretch_dx(d_stretch_u_dx, d_stretch_v_dx, ref_plane_coordinate);

    Eigen::Vector9d triangle_dot;
    triangle_dot <<1, 1, 1, 1, 1, 1, 1, 1, 1;

    Eigen::Vector9d shear_force;
    Eigen::Matrix99d d_shear_force_dx;
    Eigen::Vector9d shear_damping_force ;
    Eigen::Matrix99d d_shear_damping_force_dx;
    Eigen::Matrix99d d_shear_damping_force_dv;
    compute_shear_force(shear_force, d_shear_force_dx, shear_damping_force, d_shear_damping_force_dx, d_shear_damping_force_dv, 
                        triangle_dot, stretch, d_stretch_u_dx, d_stretch_v_dx,
                        1, 1, area);

    Eigen::Vector9d expected_shear_force;
    expected_shear_force << 1, 1, 0, 3, -1, 0, -4, 0, 0; 
    Eigen::Matrix99d expected_d_shear_force_dx; 
    expected_d_shear_force_dx << 
    -0.75,  0.25,    0,  0.75, -0.25,    0,    0,    0,    0,
     0.25, -0.75,    0,  0.75, -0.25,    0, -1.0,  1.0,    0,
        0,     0, -1.0,     0,     0,    0,    0,    0,  1.0,
     0.75,  0.75,    0,  3.25, -0.75,    0, -4.0,    0,    0,
    -0.25, -0.25,    0, -0.75,  1.25,    0,  1.0, -1.0,    0,
        0,     0,    0,     0,     0,  1.0,    0,    0, -1.0,
        0,  -1.0,    0,  -4.0,   1.0,    0,  4.0,    0,    0,
        0,   1.0,    0,     0,  -1.0,    0,    0,    0,    0,
        0,     0,  1.0,     0,     0, -1.0,    0,    0,    0;
    assertTrue(expected_shear_force == shear_force);
    assertTrue(expected_d_shear_force_dx == d_shear_force_dx);
}

void test_bend_force() {
    std::cout << "TEST test_bend_force \n";

    Eigen::Vector9d triangle_vertices1;
    Eigen::Vector9d triangle_vertices2;
    triangle_vertices1 << 1, 1, 3, 4, 0, 3, 2, 2, 3;
    triangle_vertices2 << 4, 0, 3, 1, 1, 3, 2, 0, 1;
    Eigen::Vector9d triangle_dot1;
    Eigen::Vector9d triangle_dot2;
    triangle_dot1 << 0, 0, 0, 0, 0, 0, 0, 0, 0;
    triangle_dot2 << 0, 0, 0, 0, 0, 0, 0, 0, 0;
  
    Eigen::Vector12d bend_force;
    Eigen::Matrix1212d d_bend_force_dx;
    Eigen::Vector12d bend_damping_force;
    Eigen::Matrix1212d d_bend_damping_force_dx;
    Eigen::Matrix1212d d_bend_damping_force_dv;
    double output;
    compute_bend_force(bend_force, d_bend_force_dx, bend_damping_force,
                       d_bend_damping_force_dx, d_bend_damping_force_dv,
                       triangle_dot1, triangle_dot2, triangle_vertices1, triangle_vertices2, 
                        1, 1);
    Eigen::Vector12d expected;
    expected << 
        -0.21811418457587243824562153645682,
        -0.29081891276782991766082871527576,
        0.90880910239946849269008973523675,
        0.03635236409597873970760358940947,
        -0.25446654867185117795322512586629,
        0.27264273071984054780702692057103,
                                        0,
                                        0,
        -0.99969001263941534195909870876043,
        0.18176182047989369853801794704735,
        0.54528546143968109561405384114205,
        -0.18176182047989369853801794704735;

    Eigen::Matrix1212d expected_dx;
    expected_dx << 
 -0.13548595250238270115026424938189, -0.059473389683247802508340367810951,    0.59481843883991020454910403026576,  0.071050810878368769801848827442942,   0.04792978529109970700111204904146,  -0.22143047340379307541890827442444,                                   0,                                   0,  -0.30895282381210319778178033390237,  0.064435141624013931348415421938947,  0.011543604392148095507228318769491,  -0.064435141624013931348415421938947,
-0.059473389683247802508340367810951,    0.2842257880487903270649154036801,    0.33868670058681269305376050606931, -0.026440132482104105956213528107645, 0.0033191068948350431554767497061634,  0.083429828128054438062326357116057,                                   0,                                   0,  -0.33620300654951522265153296726677,  0.085913522165351908464553895918596,  -0.28754489494362537022039215338626,  -0.085913522165351908464553895918596,
  0.59481843883991020454910403026576,   0.33868670058681269305376050606931,     0.4570432389173240457677230520407,  -0.09913640647331836742485067171096,  -0.03309243660634622217811708621268,   0.11530155321760996990575476196653, -0.49984500631970767097954935438021, -0.49984500631970767097954935438021,  -0.56818181818181818181818181818182, 0.0041629739531158338552959958254121,   0.19425074233924120010390593452359, -0.0041629739531158338552959958254121,
 0.071050810878368769801848827442942, -0.026440132482104105956213528107645,   -0.09913640647331836742485067171096,  0.030569289632247068025229383070558,   0.14953861353405792089943021260079, -0.084269468085963619788860585627493,                                   0,                                   0,   0.17266668428861299865564202034863,   -0.1016201005106158378270782105135,  -0.12309848105195381494321668449315,   0.010739190270668988558069236989825,
  0.04792978529109970700111204904146, 0.0033191068948350431554767497061634,   -0.03309243660634622217811708621268,   0.14953861353405792089943021260079,  0.043800628140956744932096194078547,   0.18092218052198451681148371853591,                                   0,                                   0, -0.072655412020955374726881973394463,  -0.19746839882515762790054226164225,  -0.04711973503579178808757294378471,  -0.075174331894682919906484658928772,
 -0.22143047340379307541890827442444,  0.083429828128054438062326357116057,    0.11530155321760996990575476196653, -0.084269468085963619788860585627493,   0.18092218052198451681148371853591,  0.020049520326891495088684992826171,  0.24992250315985383548977467719011, -0.24992250315985383548977467719011,  -0.17045454545454545454545454545455,  0.055777438329902859717994182861829, -0.014429505490185119384035398461864,   0.035103471910043989551014790661846,
                                   0,                                    0,   -0.49984500631970767097954935438021,                                    0,                                    0,   0.24992250315985383548977467719011,                                   0,                                   0,   0.24992250315985383548977467719011,                                    0,                                    0,                                     0,
                                   0,                                    0,   -0.49984500631970767097954935438021,                                    0,                                    0,  -0.24992250315985383548977467719011,                                   0,                                   0,   0.74976750947956150646932403157032,                                    0,                                    0,                                     0,
 -0.30895282381210319778178033390237,  -0.33620300654951522265153296726677,   -0.56818181818181818181818181818182,   0.17266668428861299865564202034863, -0.072655412020955374726881973394463,  -0.17045454545454545454545454545455,  0.24992250315985383548977467719011,  0.74976750947956150646932403157032,                                0.625,  -0.11363636363636363636363636363636,  -0.34090909090909090909090909090909,    0.11363636363636363636363636363636,
 0.064435141624013931348415421938947,  0.085913522165351908464553895918596,  0.0041629739531158338552959958254121,   -0.1016201005106158378270782105135,  -0.19746839882515762790054226164225,  0.055777438329902859717994182861829,                                   0,                                   0,  -0.11363636363636363636363636363636,  0.037184958886601906478662788574553,   0.11155487665980571943598836572366,   0.053695951353344942790346184949123,
 0.011543604392148095507228318769491,  -0.28754489494362537022039215338626,    0.19425074233924120010390593452359,  -0.12309848105195381494321668449315,  -0.04711973503579178808757294378471, -0.014429505490185119384035398461864,                                   0,                                   0,  -0.34090909090909090909090909090909,   0.11155487665980571943598836572366,   0.33466462997941715830796509717097,    0.16108785406003482837103855484737,
-0.064435141624013931348415421938947, -0.085913522165351908464553895918596, -0.0041629739531158338552959958254121,  0.010739190270668988558069236989825, -0.075174331894682919906484658928772,  0.035103471910043989551014790661846,                                   0,                                   0,   0.11363636363636363636363636363636,  0.053695951353344942790346184949123,   0.16108785406003482837103855484737,    -0.1445768615932917920593551584728;
    bool same = true;
    for(int i = 0; i < 12; i++) {
        if (abs(expected(i) - bend_force(i)) > 1e-15) {
            // allow small rounding error
            same = false;
            break;
        }
        for (int j = 0; j< 12; j++) {
            if (abs(expected_dx(i, j) - d_bend_force_dx(i, j)) > 1e-15) {
                // allow small rounding error
                same = false;
                break;
            }
        }
    }
    assertTrue(same);
}