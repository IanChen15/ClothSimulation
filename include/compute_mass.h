#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <EigenTypes.h>
//Input:
//  F - the mx3 triangle connectivity matrix. Each row contains to indices into V that indicate a spring between those vertices.
//  a0 - the mx1 vector of undeformed triangle areas
//  density - density of the cloth
//Output:
//  M - The 3n * 3n diaganol matrix containing the mass of individual vertex
void compute_mass(Eigen::SparseMatrixd &M,
                  Eigen::Ref<const Eigen::MatrixXi> F,
                  Eigen::Ref<const Eigen::VectorXd> a0,
                  double density);