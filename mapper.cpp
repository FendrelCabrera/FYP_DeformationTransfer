#include <vector>
#include <algorithm>
#include <utility>
#include <map>
#include <igl/heat_geodesics.h>

class Mapper
{
private:
    std::vector<std::pair<int, int>> fmap, rmap;
    std::vector<int> cmap;
    Eigen::MatrixXd &V1, &V2, dm1, dm2;
    Eigen::MatrixXi &F1, &F2;
    igl::HeatGeodesicsData<double> heat_data1, heat_data2;
    bool latestCmap = true;

    void addCol(int mid, int vid)
    {
        Eigen::VectorXd d;
        // Eigen::VectorXi VS, FS, VT, FT;
        // FS.resize(1);
        // FS << key;

        Eigen::VectorXi VS;
        VS.resize(1);

        if (mid == 0)
        {
            // FT.setLinSpaced(F1.rows(), 0, F1.rows() - 1);
            // igl::exact_geodesic(V1, F1, VS, FS, VT, FT, d);

            igl::heat_geodesics_solve(heat_data1, (Eigen::VectorXi(1, 1) << vid).finished(), d);
            std::cout << "DV1 size: " << d.rows() << " * " << d.cols() << std::endl;

            dm1.conservativeResize(V1.rows(), dm1.cols() + 1);
            dm1.col(dm1.cols() - 1) = d;
        }
        else
        {
            // FT.setLinSpaced(F2.rows(), 0, F2.rows() - 1);
            // igl::exact_geodesic(V2, F2, VS, FS, VT, FT, d);

            igl::heat_geodesics_solve(heat_data2, (Eigen::VectorXi(1, 1) << vid).finished(), d);
            std::cout << "DV2 size: " << d.rows() << " * " << d.cols() << std::endl;

            dm2.conservativeResize(V2.rows(), dm2.cols() + 1);
            dm2.col(dm2.cols() - 1) = d;
        }

        latestCmap = false;
    }

    void delCol(int mid, int index)
    {
        if (mid == 0)
        {
            dm1.block(0, index, dm1.rows(), dm1.cols() - index - 1) = dm1.rightCols(dm1.cols() - index - 1);
            dm1.conservativeResize(dm1.rows(), dm1.cols() - 1);
        }
        else
        {
            dm2.block(0, index, dm2.rows(), dm2.cols() - index - 1) = dm2.rightCols(dm2.cols() - index - 1);
            dm2.conservativeResize(dm2.rows(), dm2.cols() - 1);
        }
        latestCmap = false;
    }

    void updateCM()
    {
        if (fmap.size() == 0)
            return;

        std::vector<std::vector<std::pair<double, int>>>
            g1(fmap.size(), std::vector<std::pair<double, int>>(0)),
            g2(rmap.size(), std::vector<std::pair<double, int>>(0));

        if (cmap.size() != V1.rows())
        {
            cmap.resize(V1.rows(), -1);
        }

        Eigen::Index minIndex;
        double dist;

        for (int i = 0; i < V1.rows(); i++)
        {
            // find minimum value in the ith row of dm1 and the index at which it occurs in that row
            dist = dm1.row(i).minCoeff(&minIndex);
            g1[minIndex].push_back({dist, i});
            cmap[i] = minIndex;
        }

        for (auto &x : g1)
            std::sort(x.begin(), x.end());

        for (int i = 0; i < V2.rows(); i++)
        {
            dist = dm2.row(i).minCoeff(&minIndex);
            g2[minIndex].push_back({dist, i});
        }

        for (auto &x : g2)
            std::sort(x.begin(), x.end());

        for (int i = 0; i < V1.rows(); i++)
        {
            auto [closestMarked, coMarked] = fmap[cmap[i]];

            auto i1 = std::find_if(g1[cmap[i]].begin(), g1[cmap[i]].end(),
                                   [i](const std::pair<double, int> &p)
                                   { return p.second == i; });
            auto percent = (i1 - g1[cmap[i]].begin()) / (double)g1[cmap[i]].size();

            auto coIndex = std::find_if(rmap.begin(), rmap.end(),
                                        [coMarked](const std::pair<int, int> &p)
                                        { return p.first == coMarked; }) -
                           rmap.begin();

            auto coFace = g2[coIndex][percent * g2[coIndex].size()].second;
            cmap[i] = coFace;
        }
    };

public:
    Mapper(Eigen::MatrixXd &V1, Eigen::MatrixXi &F1, Eigen::MatrixXd &V2, Eigen::MatrixXi &F2) : V1(V1), V2(V2), F1(F1), F2(F2)
    {
        igl::heat_geodesics_precompute(V1, F1, heat_data1);
        std::cout << "M1 precompute complete" << std::endl;
        igl::heat_geodesics_precompute(V2, F2, heat_data2);
        std::cout << "M2 precompute complete" << std::endl;
    }

    int findClosestVertex(
        int mid,
        int fid, Eigen::Vector3f bc)
    {
        Eigen::MatrixXd V;
        Eigen::MatrixXi F;
        if (mid == 0)
        {
            V = V1;
            F = F1;
        }
        else
        {
            V = V2;
            F = F2;
        }

        Eigen::RowVector3d c = bc(0) * V.row(F(fid, 0)) +
                               bc(1) * V.row(F(fid, 1)) +
                               bc(2) * V.row(F(fid, 2));
        int cid = 0;
        Eigen::Vector3d(
            (V.row(F(fid, 0)) - c).squaredNorm(),
            (V.row(F(fid, 1)) - c).squaredNorm(),
            (V.row(F(fid, 2)) - c).squaredNorm())
            .minCoeff(&cid);
        return F(fid, cid);
    }

    void set(int key, int value)
    {
        auto it = std::find_if(fmap.begin(), fmap.end(),
                               [key](const std::pair<int, int> &p)
                               { return p.first == key; });

        if (it != fmap.end())
        {
            auto i2 = std::find_if(rmap.begin(), rmap.end(),
                                   [it](const std::pair<int, int> &p)
                                   { return p.first == it->second; });
            if (i2->second == 1)
            {
                delCol(1, i2 - rmap.begin());
                rmap.erase(i2);
            }
            else
                i2->second--;
            it->second = value;
        }
        else
        {
            fmap.push_back({key, value});
            addCol(0, key);
        }

        it = std::find_if(rmap.begin(), rmap.end(),
                          [value](const std::pair<int, int> &p)
                          { return p.first == value; });
        if (it != rmap.end())
        {
            it->second++;
        }
        else
        {
            rmap.push_back({value, 1});
            addCol(1, value);
        }

        // updateCM();
    }

    int get(int key)
    {
        auto it = std::find_if(fmap.begin(), fmap.end(),
                               [key](const std::pair<int, int> &p)
                               { return p.first == key; });
        return it != fmap.end() ? it->second : -1;
    }

    void erase(int key)
    {
        auto i1 = std::find_if(fmap.begin(), fmap.end(),
                               [key](const std::pair<int, int> &p)
                               { return p.first == key; });
        if (i1 != fmap.end())
        {
            auto i2 = std::find_if(rmap.begin(), rmap.end(),
                                   [i1](const std::pair<int, int> &p)
                                   { return p.first == i1->second; });
            if (i2->second == 1)
            {
                delCol(1, i2 - rmap.begin());
                rmap.erase(i2);
            }
            else
                i2->second--;
            delCol(0, i1 - fmap.begin());
            fmap.erase(i1);
        }

        // updateCM();
    }

    void clear()
    {
        fmap.clear();
        rmap.clear();
        dm1.resize(V1.rows(), 0);
        dm2.resize(V2.rows(), 0);
        std::fill(cmap.begin(), cmap.end(), -1);
        latestCmap = true;
    }

    int vcount(int value)
    {
        auto it = std::find_if(rmap.begin(), rmap.end(),
                               [value](const std::pair<int, int> &p)
                               { return p.first == value; });
        return it != rmap.end() ? it->second : 0;
    }

    int corr(int key)
    {
        if (!latestCmap)
            updateCM();
        if (cmap.size() != V1.rows())
            cmap.resize(V1.rows(), -1);
        return cmap[key];
    }

    Eigen::MatrixXd getFMarked()
    {
        Eigen::MatrixXd result(fmap.size(), 3);
        for (int i = 0; i < fmap.size(); i++)
        {
            result.row(i) = V1.row(fmap[i].first);
        }
        return result;
    }

    Eigen::MatrixXd getRMarked()
    {
        Eigen::MatrixXd result(rmap.size(), 3);
        for (int i = 0; i < rmap.size(); i++)
        {
            result.row(i) = V2.row(rmap[i].first);
        }
        return result;
    }
};
