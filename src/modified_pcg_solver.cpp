#include <modified_pcg_solver.h>
#include <iostream>

void filter(Eigen::VectorXd &a, const std::vector<unsigned int> &fixed_points);
void precondition_matrix(Eigen::SparseMatrixd &P, Eigen::SparseMatrixd &P_inv, const Eigen::SparseMatrixd &A);
void modified_pcg_solver(Eigen::VectorXd &d_qdot, 
                         const Eigen::SparseMatrixd &A, 
                         Eigen::Ref<const Eigen::VectorXd> b,
                         const std::vector<unsigned int> &fixed_points) {
    d_qdot.resize(b.rows());
    Eigen::VectorXd z;
    z.resize(d_qdot.rows());
    z.setZero();
    Eigen::SparseMatrixd P;
    Eigen::SparseMatrixd P_inv;
    precondition_matrix(P, P_inv, A);
    d_qdot = z;

    Eigen::VectorXd filter_b = b;
    filter(filter_b, fixed_points);
    double delta0 = filter_b.transpose() * P * filter_b;

    Eigen::VectorXd r = b - A * d_qdot;
    filter(r, fixed_points);
    
    Eigen::VectorXd c = P_inv * r;
    filter(c, fixed_points);

    double delta_new = r.transpose() * c;
    double tol = 0.0001;
    int max_it = 50;
    int i;
	for (i = 0; i < max_it && delta_new > tol * tol * delta0; i++) {
        Eigen::VectorXd q = A * c;
        filter(q, fixed_points);
        double alpha = delta_new / (c.transpose() * q);
        d_qdot = d_qdot + alpha * c;
        r = r - alpha * q;
        Eigen::VectorXd s = P_inv * r;
        double delta_old = delta_new;
        delta_new = r.transpose() * s;
        c = s + delta_new / delta_old * c;
        filter(c, fixed_points);
    } 
    // remove small values
    for(int i = 0; i < d_qdot.rows(); i++) {
        if (std::abs(d_qdot(i)) < 1e-16) {
            d_qdot(i) = 0;
        }
    }
}

void filter(Eigen::VectorXd &a, const std::vector<unsigned int> &fixed_points) {
    for (int i = 0; i < fixed_points.size(); i ++) {
        unsigned int fixed_point = fixed_points[i];
        a(3 * fixed_point) = 0;
        a(3 * fixed_point + 1) = 0;
        a(3 * fixed_point + 2) = 0;
    }
}

void precondition_matrix(Eigen::SparseMatrixd &P, Eigen::SparseMatrixd &P_inv, const Eigen::SparseMatrixd &A) {
    P.resize(A.rows(), A.rows());
    P_inv.resize(A.rows(), A.rows());
    P.setZero();
    P_inv.setZero();
    std::vector<Eigen::Triplet<double>> update_list;  
    std::vector<Eigen::Triplet<double>> update_list_inv;  
    for(int i = 0; i < A.rows(); i++) {
        update_list.push_back(Eigen::Triplet<double>(i, i, 1.0/A.coeff(i, i)));  
        update_list_inv.push_back(Eigen::Triplet<double>(i, i, A.coeff(i, i)));  
    }
    P.setFromTriplets(update_list.begin(), update_list.end());  
    P_inv.setFromTriplets(update_list_inv.begin(), update_list_inv.end());  
}