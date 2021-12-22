#include <dV_cloth_gravity_dq.h>

void dV_cloth_gravity_dq(Eigen::VectorXd &fg, Eigen::SparseMatrixd &M, Eigen::Ref<const Eigen::Vector3d> g) {
    // Expand the gravity for every point
    Eigen::VectorXd G = g.replicate(M.rows()/3, 1);
    fg = -M * G;
}
