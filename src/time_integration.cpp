#include <time_integration.h>
#include <iostream>

void time_integration(Eigen::VectorXd &q, 
                      Eigen::VectorXd &qdot, 
                      const Eigen::SparseMatrixd &M, 
                      Eigen::Ref< const Eigen::VectorXd> f, 
                      const Eigen::SparseMatrixd &df_dx, 
                      const Eigen::SparseMatrixd &df_dv, 
                      const std::vector<unsigned int> &fixed_points,
                      double h) {
    Eigen::SparseMatrixd A = M - h * df_dv - h * h * df_dx;
    Eigen::VectorXd b = h * (f + h * df_dx * qdot);
    Eigen::VectorXd d_qdot;
    modified_pcg_solver(d_qdot, A, b, fixed_points);

    qdot = qdot + d_qdot;
    q = q + h * qdot;
};

