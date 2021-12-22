#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <EigenTypes.h>

#include <modified_pcg_solver.h>

void time_integration(Eigen::VectorXd &q, 
                      Eigen::VectorXd &qdot, 
                      const Eigen::SparseMatrixd &M, 
                      Eigen::Ref< const Eigen::VectorXd> f, 
                      const Eigen::SparseMatrixd &df_dx, 
                      const Eigen::SparseMatrixd &df_dv,
                      const std::vector<unsigned int> &fixed_points,
                      double h);