#ifndef ASSIGNMENT_SETUP_H
#define ASSIGNMENT_SETUP_H

//Assignment 4 code 
#include <igl/readMESH.h>
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/readOFF.h>
#include <read_tetgen.h>
#include <igl/boundary_facets.h>
#include <igl/volume.h>

//assignment files for implementing simulation and interaction
#include <visualization.h>
#include <init_state.h>
#include <find_max_vertices.h>

#include <dV_cloth_gravity_dq.h>
#include <dV_spring_particle_particle_dq.h>
#include <compute_mass.h>
#include <assemble_forces.h>
#include <time_integration.h>

extern double K_BEND;
extern double K_STRETCH;
extern double K_SHEAR;
extern double K_DAMP;
extern double bu;
extern double bv;
#define GRAVITY 9.81
#define DENSITY 0.1

//Variable for geometry
Eigen::MatrixXd V, V_skin; //vertices of simulation mesh //this will hold all individual pieces of cloth, I'll load some offsets
Eigen::MatrixXi F, F_skin; //faces of simulation mesh

Eigen::MatrixXd V_sphere, V_sphere_skin; //vertices of simulation mesh //this will hold all individual pieces of cloth, I'll load some offsets
Eigen::MatrixXi F_sphere, F_sphere_skin; //faces of simulation mesh

Eigen::SparseMatrixd N;

//material parameters
double density = DENSITY;

//BC
std::vector<unsigned int> fixed_point_indices;

//mass matrix
Eigen::SparseMatrixd M; //mass matrix
Eigen::VectorXd a0; //areas
Eigen::MatrixXd uv_coords; // n * 6, which are 2d coordinate of triangle vertices
std::vector<std::vector<int>> adjacency_list;

Eigen::VectorXd gravity;

std::vector<std::pair<Eigen::Vector3d, unsigned int>> spring_points; //need this for interaction 

//collision detection stuff
bool collision_detection_on = false; 
std::vector<unsigned int> collision_indices;
std::vector<Eigen::Vector3d> collision_normals;

//selection spring
double k_selected = 1e5;
inline void simulate(Eigen::VectorXd &q, Eigen::VectorXd &qdot, double dt, double t) {  

    double V_ele, T_ele, KE,PE;

    spring_points.clear();

    //Interaction spring
    Eigen::Vector3d mouse;
    Eigen::Vector6d dV_mouse;
    double k_selected_now = (Visualize::is_mouse_dragging() ? k_selected : 0.);
    
    for(unsigned int pickedi = 0; pickedi < Visualize::picked_vertices().size(); pickedi++) {   
        spring_points.push_back(std::make_pair((q).segment<3>(3*Visualize::picked_vertices()[pickedi]) + Visualize::mouse_drag_world() + Eigen::Vector3d::Constant(1e-6),3*Visualize::picked_vertices()[pickedi]));
    }

    Eigen::VectorXd f;
    Eigen::SparseMatrixd df_dx;
    Eigen::SparseMatrixd df_dv;

    assemble_forces(f, df_dx, df_dv, q, qdot, F, uv_coords,  adjacency_list, a0);
     f -= gravity;

    for(unsigned int pickedi = 0; pickedi < spring_points.size(); pickedi++) {
        dV_spring_particle_particle_dq(dV_mouse, spring_points[pickedi].first, (q).segment<3>(spring_points[pickedi].second), 0.0, k_selected_now);
        f.segment<3>(3*Visualize::picked_vertices()[pickedi]) -= dV_mouse.segment<3>(3) * 0.0004;
        std::cout <<dV_mouse.segment<3>(3).transpose() * 0.0004 << " " << 3*Visualize::picked_vertices()[pickedi] << "\n";
    }
    
    time_integration(q, qdot, M, f, df_dx, df_dv, fixed_point_indices, dt);// {    // KE = PE = 0.;
 }

inline void draw(Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot, double t) {

    //update vertex positions using simulation
    Visualize::update_vertex_positions(0, q);

}

bool key_down_callback(igl::opengl::glfw::Viewer &viewer, unsigned char key, int modifiers) {
    if(key == 'B') {
        if (K_BEND == 0.001) {
            K_BEND = 0.01;   
        } else if (K_BEND == 0.01) {
            K_BEND = 0.05;
        } else {
            K_BEND = 0.001;
        }
        std::cout << "K_BEND updated to: " << K_BEND << "\n";
    }
    if(key == 'S') {
        if (K_STRETCH == 5000.0) {
            K_STRETCH = 2000.0;   
        } else if (K_STRETCH == 2000.0) {
            K_STRETCH = 1000.0;
        } else {
            K_STRETCH = 5000.0;
        }
        std::cout << "K_STRETCH updated to: " << K_STRETCH << "\n";
    }
    if(key == 'R') {
        if (K_SHEAR == 5000.0) {
            K_SHEAR = 500.0;   
        } else if (K_SHEAR == 500.0) {
            K_SHEAR = 300.0;
        } else {
            K_SHEAR = 5000.0;
        }
        std::cout << "K_SHEAR updated to: " << K_SHEAR << "\n";
    }
    if(key == 'D') {
        if (K_DAMP == 0.2) {
            K_DAMP = 0.4;   
        } else if (K_DAMP == 0.4) {
            K_DAMP = 0.6;
        } else {
            K_DAMP = 0.2;
        }
        std::cout << "K_DAMP updated to: " << K_DAMP << "\n";
    }
    return false;
}
inline void assignment_setup(int argc, char **argv, Eigen::VectorXd &q, Eigen::VectorXd &qdot) {

    //load geometric data 
    igl::readOBJ("../data/square_cloth.obj", V, F);
    igl::readOBJ("../data/sphere.obj", V_sphere, F_sphere);

    //setup simulation 
    init_state(q,qdot,V);

    //add geometry to scene
    V_skin = V;
    F_skin = F;
    N.resize(V.rows(), V.rows());
    N.setIdentity();

    Visualize::add_object_to_scene(V,F, V_skin, F_skin, N, Eigen::RowVector3d(244,165,130)/255.);

    //add collision sphere to scene 
    V_sphere_skin = V_sphere;
    F_sphere_skin = F_sphere;
    N.resize(V.rows(), V.rows());
    N.setIdentity();

    Visualize::add_object_to_scene(V_sphere,F_sphere, V_sphere_skin, F_sphere_skin, N, Eigen::RowVector3d(244,165,130)/255.);
    Visualize::set_visible(1, collision_detection_on);

    //compute uv_coords and area of each face;
    uv_coords.resize(F.rows(), 6);
    a0.resize(F.rows(),1);
    for(unsigned int ii=0; ii<F.rows(); ++ii) {
        Eigen::Vector6d uv_tmp;
        uv_tmp(0) = V(F(ii, 0), 0);
        uv_tmp(1) = V(F(ii, 0), 1);
        uv_tmp(2) = V(F(ii, 1), 0);
        uv_tmp(3) = V(F(ii, 1), 1);
        uv_tmp(4) = V(F(ii, 2), 0);
        uv_tmp(5) = V(F(ii, 2), 1);
        uv_coords.row(ii) = uv_tmp;
        //Heron's formula for area
        double side_a = (V.row(F(ii,1)) - V.row(F(ii,0))).norm();
        double side_b = (V.row(F(ii,2)) - V.row(F(ii,1))).norm();
        double side_c = (V.row(F(ii,0)) - V.row(F(ii,2))).norm();
        double s = (side_a+side_b+side_c)/2.;
        a0(ii) = sqrt(s*(s-side_a)*(s-side_b)*(s-side_c));
     }

    // Compute adjacency triangles
    adjacency_list.resize(F.rows());
    for (int i = 0; i < F.rows(); i++) {
        std::vector<int> temp1 =  {F(i, 0), F(i, 1), F(i, 2)};

        for(int j = i + 1; j < F.rows(); j++) {
            std::vector<int> temp2 =  {F(j, 0), F(j, 1), F(j, 2)};
            int count = 0;
            for (int num = 0; num < 3; num++) {
                for(int num2 = 0; num2 < 3; num2++) {
                    if (temp1.at(num) == temp2.at(num2)) {
                        count++;
                    }
                }
            }
            if (count >= 2) {
                adjacency_list.at(i).push_back(j);
            }
        }
    }

    //Mass Matrix
    M.resize(q.rows(), q.rows());
    compute_mass(M, F, a0, density);
    
    //constant gravity vector
    gravity.resize(q.rows(),1);
    dV_cloth_gravity_dq(gravity, M, Eigen::Vector3d(0,-GRAVITY, 0));

       //should be max verts for cloth simulation
    // Two configuration
    if (argc > 1 && strcmp(argv[1], "1") == 0) {
        // this will have two corners of cloths fixed
        int max = 0;
        int min = 0;
        for (int i = 0; i < V.rows(); i++) {
            if (V(i, 0) > V(max, 0) && V(i, 1) > V(max, 1)) {
                max = i;
            }
            if (V(i, 0) < V(min, 0) && V(i, 1) < V(min, 1)) {
                min = i;
            }
        }
        fixed_point_indices.push_back(min);
        fixed_point_indices.push_back(max);
            
        for (int i = 1; i < q.size(); i+= 3) {
            q(i+1) = q(i);
            q(i)  = 0;
        }
    } else {
        // this will have top line of the cloth fixed
        find_max_vertices(fixed_point_indices, V, 0.001);
    }
    std::cout << "  D       Toggle Damp coefficient \n";
    std::cout << "  B       Toggle Bend coefficient \n";
    std::cout << "  S       Toggle Stretch coefficient \n";
    std::cout << "  R       Toggle Shear coefficient \n";
 
    Visualize::viewer().callback_key_down = key_down_callback;

}

#endif

