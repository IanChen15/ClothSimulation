#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <EigenTypes.h>
void modified_pcg_solver(Eigen::VectorXd &d_qdot, 
                         const Eigen::SparseMatrixd &A, 
                         Eigen::Ref<const Eigen::VectorXd> b,
                         const std::vector<unsigned int> &fixed_points);