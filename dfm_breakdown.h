// #include <Eigen/Dense>

Eigen::Vector3d estimate_force_vector(
    const Eigen::MatrixXd &V_init,
    const Eigen::MatrixXd &V_final,
    const Eigen::RowVector3d &x0,
    double epsilon,
    double mu = 1.0,
    double nu = 0.4);

std::vector<std::pair<int, Eigen::Vector3d>> find_impacts_and_forces(
    const Eigen::MatrixXd &V_init,
    const Eigen::MatrixXd &V_final,
    double epsilon = 1e-1,
    double mu = 1.0,
    double nu = 0.4);

Eigen::MatrixXd impactPoints(
    Eigen::MatrixXd Ov,
    Eigen::MatrixXd V);

void printError(
    Eigen::MatrixXd V1,
    Eigen::MatrixXd V2);