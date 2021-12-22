#include <dV_spring_particle_particle_dq.h>

void dV_spring_particle_particle_dq(Eigen::Ref<Eigen::Vector6d> f, Eigen::Ref<const Eigen::Vector3d> q0,  Eigen::Ref<const Eigen::Vector3d>     q1, double l0, double stiffness) {
    Eigen::Vector6d q;  
    q << q0, q1;  
      
    Eigen::Matrix36d B;  
    B << -1, 0, 0, 1, 0, 0,  
         0, -1, 0, 0, 1, 0,  
         0, 0, -1, 0, 0, 1;  
    Eigen::Vector3d diff = (q0 - q1);  
    // This is the same as q^t * B^t * B * q  
    double distance = diff.dot(diff);  
  
    f = stiffness * (1.0 - l0 / sqrt(distance)) * (B.transpose() * B * q);  
}