#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <EigenTypes.h>

/**
 * Precondition: Element1 and element2 are adjacent traingle.
 * Order the elements in the following order
 *       t02 
 *      /   \
 *   t00 - - t01
 *   t11 - - t10
 *      \   /
 *       t12
 * 
 * Input: 
 *   - element1, element2
 * Outut: 
 *   - ordered_element1, ordered_element2 
 **/
void order_triangles(Eigen::Vector3i &ordered_element1, 
                     Eigen::Vector3i &ordered_element2,
                     Eigen::Ref<const Eigen::Vector3i> element1,
                     Eigen::Ref<const Eigen::Vector3i> element2);

/**Input:
* 
*  q - generalized coordinates for the FEM system
*  qdot - generalized velocity for the FEM system
*  F - the mx3 triangle connectivity matrix. Each row contains to indices into V that indicate a spring between those vertices.
*  uv_coords - the nx6 matrix of undeformed vertex positions. Each row is a single undeformed vertex position.
*  adjacency_list - the n size vector which contains a list of index to the adjecent triangle. It will only contains 
*                    triangles with index after the current triangle
*  a0 - the mx1 vector of undeformed triangle areas
*Output:
*  f - the vector 3xn vector of forces acting on each node of the mass-spring system
*  df_dx - the vector 3n x 3n gradient matrix of f with respect to x
*  df_dv - the vector 3n x 3n gradient matrix of f with respect to v
**/
void assemble_forces(Eigen::VectorXd &f, 
                     Eigen::SparseMatrixd &df_dx, 
                     Eigen::SparseMatrixd &df_dv,
                     Eigen::Ref<const Eigen::VectorXd> q,
                     Eigen::Ref<const Eigen::VectorXd> qdot,
                     Eigen::Ref<const Eigen::MatrixXi> F,
                     Eigen::Ref<const Eigen::MatrixXd> uv_coords,
                     const std::vector<std::vector<int>> &adjacency_list,
                     Eigen::Ref<const Eigen::VectorXd> a0);