#include <compute_mass.h>
void compute_mass(Eigen::SparseMatrixd &M,
                  Eigen::Ref<const Eigen::MatrixXi> F,
                  Eigen::Ref<const Eigen::VectorXd> a0,
                  double density) { 
    M.setZero();  
    std::vector<Eigen::Triplet<double>> update_list;  
    Eigen::VectorXd masses;
    masses.resize(M.rows()/3);
    masses.setZero() ;
    for (int triangle_index = 0; triangle_index < F.rows(); triangle_index++) {  
        Eigen::Vector3i element = F.row(triangle_index);
        double mass = a0(triangle_index) * density / 3.0;
        masses(element(0)) += mass;
        masses(element(1)) += mass;
        masses(element(2)) += mass;
    }

    for (int i = 0; i < masses.rows(); i++) {
        update_list.push_back(Eigen::Triplet<double>(3*i, 3*i, masses(i)));  
        update_list.push_back(Eigen::Triplet<double>(3*i+1, 3*i+1, masses(i)));  
        update_list.push_back(Eigen::Triplet<double>(3*i+2, 3*i+2, masses(i)));  
    }

    M.setFromTriplets(update_list.begin(), update_list.end());  
}