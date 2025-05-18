#include <Eigen/Dense>
#include <Eigen/Core>
#include <vector>
#include <utility>
#include <iostream>
#include <cmath>
#include <set>

// Check if a vertex is a local maximum in displacement magnitude
bool is_local_maximum(const Eigen::VectorXd &magnitudes,
                      const std::vector<std::vector<int>> &adjacency,
                      int idx)
{
    for (int neighbor : adjacency[idx])
    {
        if (magnitudes[neighbor] > magnitudes[idx])
            return false;
    }
    return true;
}

Eigen::Vector3d estimate_force_vector(
    const Eigen::MatrixXd &V_init,
    const Eigen::MatrixXd &V_final,
    const Eigen::RowVector3d &x0,
    double epsilon,
    double mu,
    double nu)
{
    const int n = V_init.rows();
    Eigen::MatrixXd A(3 * n, 3);
    Eigen::VectorXd b(3 * n);

    auto phi = [&](const Eigen::RowVector3d &x)
    {
        Eigen::RowVector3d r = x - x0;
        double r_norm = std::sqrt(r.squaredNorm() + epsilon * epsilon);
        return (1.0 / (8.0 * M_PI * mu * (1.0 - nu))) *
               ((3.0 - 4.0 * nu) / r_norm + (epsilon * epsilon) / (r_norm * r_norm * r_norm));
    };

    for (int i = 0; i < n; ++i)
    {
        double scalar = phi(V_init.row(i));
        Eigen::Matrix3d Ki = scalar * Eigen::Matrix3d::Identity();
        A.block<3, 3>(3 * i, 0) = Ki;
        b.segment<3>(3 * i) = (V_final.row(i) - V_init.row(i)).transpose();
    }

    return A.colPivHouseholderQr().solve(b); // Least-squares force vector
}

std::vector<std::vector<int>> build_vertex_adjacency(const Eigen::MatrixXi &F, int num_vertices)
{
    std::vector<std::set<int>> adj_sets(num_vertices);
    for (int i = 0; i < F.rows(); ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            int v0 = F(i, j);
            int v1 = F(i, (j + 1) % 3);
            int v2 = F(i, (j + 2) % 3);
            adj_sets[v0].insert(v1);
            adj_sets[v0].insert(v2);
        }
    }

    std::vector<std::vector<int>> adjacency(num_vertices);
    for (int i = 0; i < num_vertices; ++i)
        adjacency[i] = std::vector<int>(adj_sets[i].begin(), adj_sets[i].end());

    return adjacency;
}

Eigen::MatrixXd impactPoints(
    Eigen::MatrixXd Ov,
    Eigen::MatrixXd V)
{
    Eigen::MatrixXd P(Ov.rows(), 3);

    // Estimate the size of the mesh
    double meshSize = (Ov.colwise().maxCoeff() - Ov.colwise().minCoeff()).norm();
    double threshold = 0.15 * meshSize; // Set threshold as a fraction of the mesh size
    std::cout << "Mesh size: " << meshSize << ", Threshold: " << threshold << std::endl;

    double diffThreshold = 0.05 * meshSize;
    int count = 0;
    std::vector<std::pair<int, double>> differences;
    for (int i = 0; i < Ov.rows(); i++)
    {
        double diff = (V.row(i) - Ov.row(i)).norm();
        if (diff < diffThreshold)
            continue;
        differences.emplace_back(i, diff);
    }
    std::cout << "Number of points with significant displacement: " << differences.size() << std::endl;

    std::sort(differences.begin(), differences.end(),
              [](const std::pair<int, double> &a, const std::pair<int, double> &b)
              {
                  return a.second > b.second;
              });

    std::vector<int> selectedIndices;
    for (const auto &pair : differences)
    {
        int idx = pair.first;
        Eigen::RowVector3d candidate = Ov.row(idx);
        bool isFarEnough = true;

        for (int selectedIdx : selectedIndices)
        {
            if ((candidate - Ov.row(selectedIdx)).norm() < threshold)
            {
                isFarEnough = false;
                break;
            }
        }

        if (isFarEnough)
        {
            selectedIndices.push_back(idx);
            P.row(count++) = candidate;
            if (count >= 5) // Limit to 10 points
                break;
        }
    }

    P.conservativeResize(count, 3);
    return P;
}

std::vector<std::pair<int, Eigen::Vector3d>> find_impacts_and_forces(
    const Eigen::MatrixXd &Ov,
    const Eigen::MatrixXd &V,
    double epsilon,
    double mu,
    double nu)
{
    // const int n = V_init.rows();
    // Eigen::VectorXd magnitudes(n);
    // for (int i = 0; i < n; ++i)
    //     magnitudes[i] = (V_final.row(i) - V_init.row(i)).norm();

    // std::vector<std::vector<int>> adjacency = build_vertex_adjacency(F, n);

    std::vector<std::pair<int, Eigen::Vector3d>> results;

    // auto P = impactPoints(V_init, V_final);
    // int n = P.rows();

    // for (int i = 0; i < n; ++i)
    // {
    //     Eigen::RowVector3d x0 = P.row(i);
    //     Eigen::Vector3d f = estimate_force_vector(V_init, V_final, x0, epsilon, mu, nu);
    //     results.emplace_back(x0, f);
    // }
    double meshSize = (Ov.colwise().maxCoeff() - Ov.colwise().minCoeff()).norm();
    double threshold = 0.15 * meshSize; // Set threshold as a fraction of the mesh size
    std::cout << "Mesh size: " << meshSize << ", Threshold: " << threshold << std::endl;

    double diffThreshold = 0.03 * meshSize;
    int count = 0;
    std::vector<std::pair<int, Eigen::Vector3d>> differences;
    for (int i = 0; i < Ov.rows(); i++)
    {
        auto diff = (V.row(i) - Ov.row(i));
        if (diff.norm() < diffThreshold)
            continue;
        differences.emplace_back(i, diff);
    }
    std::cout << "Number of points with significant displacement: " << differences.size() << std::endl;

    std::sort(differences.begin(), differences.end(),
              [](const std::pair<int, Eigen::Vector3d> &a, const std::pair<int, Eigen::Vector3d> &b)
              {
                  return a.second.norm() > b.second.norm();
              });

    std::vector<int> selectedIndices;
    for (const auto &pair : differences)
    {
        int idx = pair.first;
        Eigen::RowVector3d candidate = Ov.row(idx);
        bool isFarEnough = true;

        for (int selectedIdx : selectedIndices)
        {
            if ((candidate - Ov.row(selectedIdx)).norm() < threshold)
            {
                isFarEnough = false;
                break;
            }
        }

        if (isFarEnough)
        {
            selectedIndices.push_back(idx);
            // P.row(count++) = candidate;
            // if (count >= 5) // Limit to 10 points
            //     break;
            results.push_back({idx, pair.second});
        }
    }

    return results;
}

void printError(
    Eigen::MatrixXd V1,
    Eigen::MatrixXd V2)
{
    double sumDiff = 0, maxDiff = 0;

    for (int i = 0; i < V1.rows(); i++)
    {
        double diff = (V1.row(i) - V2.row(i)).norm();
        sumDiff += diff;
        if (diff > maxDiff)
            maxDiff = diff;
    }

    double avgDiff = sumDiff / V1.rows();
    std::cout << avgDiff << "," << maxDiff << std::endl;
}