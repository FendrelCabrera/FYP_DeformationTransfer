#define IGL_VIEWER_VIEWER_QUIET 0
#include <igl/opengl/glfw/Viewer.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/unproject.h>
#include <igl/kelvinlets.h>
#include <iostream>
#include <chrono>
#include "dfm_breakdown.h"
#include "mapper.cpp"

unsigned core1, core2, core3;
typedef enum
{
  MRK,
  CRS,
  DFM
} OperatingMode;

void applyKVL(
    Eigen::MatrixXd &V,
    Eigen::Vector3d &posStart,
    Eigen::Vector3d &forceVec,
    const igl::KelvinletParams<double> &params,
    Eigen::MatrixXd &result);

int main(int argc, char *argv[])
{
  if (argc < 3)
  {
    std::cout << "./example <sourceMesh> <targetMesh>" << std::endl;
    exit(1);
  }

  igl::opengl::glfw::Viewer viewer;
  Eigen::MatrixXd Vs1, Vs2, Vt, Ot;
  Eigen::MatrixXi Fs, Ft;

  auto brushRadius = 1.;
  auto brushType = igl::BrushType::GRAB;
  auto brush_strength = 1.;
  auto scale = 1;

  Eigen::Vector3d posStart(0, 0, 0), posEnd, projectedPoint;
  decltype(Vs1) result;

  viewer.load_mesh_from_file(argv[1]);
  viewer.load_mesh_from_file(argv[1]);
  viewer.load_mesh_from_file(argv[2]);

  Vs1 = Vs2 = viewer.data_list[0].V;
  Fs = viewer.data_list[0].F;
  Vt = Ot = viewer.data_list[2].V;
  Ft = viewer.data_list[2].F;

  std::cout << "(" << Vs1.rows() << "," << Fs.rows() << ") (" << Vt.rows() << "," << Ft.rows() << ")" << std::endl;

  OperatingMode mode = MRK;
  Mapper m(Vs1, Fs, Vt, Ft);

  brush_strength = (Vs1.colwise().maxCoeff() - Vs1.colwise().minCoeff()).norm();
  brushRadius = brush_strength / 20;
  auto bs2 = (Ot.colwise().maxCoeff() - Ot.colwise().minCoeff()).norm();
  auto br2 = bs2 / 20;

  const auto launchMenu = [&]()
  {
    std::cout << "Menu" << std::endl;
    std::cout << "1. Change brush radius" << std::endl;
    std::cout << "2. Change brush strength" << std::endl;
    std::cout << "3. Change scale" << std::endl;
    std::cout << "4. Reset deformation" << std::endl;
    std::cout << "5. Compute impact points and apply deformation" << std::endl;
    std::cout << "6. Cancel" << std::endl;

    std::cout << "Enter your choice: ";
    int choice;
    std::cin >> choice;
    switch (choice)
    {
    case 1:
      std::cout << "Brush Radius " << brushRadius << " -> ";
      std::cin >> brushRadius;
      break;
    case 2:
      std::cout << "Brush Strength " << brush_strength << " -> ";
      std::cin >> brush_strength;
      break;
    case 3:
      std::cout << "Scale " << scale << " -> ";
      std::cin >> scale;
      break;
    case 4:
      std::cout << "Resetting deformation" << std::endl;
      Vs2 = Vs1;
      viewer.data_list[1].set_vertices(Vs2);
      viewer.data_list[1].compute_normals();
      break;
    case 5:
    {
      // mode = TRF;
      // std::cout << "Click on point in right mesh" << std::endl;

      // Eigen::MatrixXd P = impactPoints(Ov, V1);
      // std::cout << "Number of vertices changed: " << P.rows() << std::endl;
      // viewer.data_list[1].add_points(P, Eigen::RowVector3d(1, 0, 0));

      // std::cout << "Impacts on " << out.size() << " vertices" << std::endl;
      // result = V1;
      // Eigen::MatrixXd t;
      // for (int i = 0; i < out.size(); i++)
      // {
      //   Eigen::Vector3d point = out[i].first.transpose();
      //   applyKVL(
      //       result,
      //       point,
      //       out[i].second,
      //       igl::KelvinletParams<double>(brushRadius, scale, brushType),
      //       t);
      //   result = t;
      // }

      // viewer.data_list[1].set_vertices(result);
      // viewer.data_list[1].compute_normals();
      // std::cout << "DT complete" << std::endl;
      // std::cout << elapsed.count() << ",";
      // printError(V1, result);

      auto start = std::chrono::high_resolution_clock::now();
      auto impactPoints = find_impacts_and_forces(
          Vs1,
          Vs2, 0.25);
      auto end = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double, std::milli> elapsed = end - start;
      std::cout << "Time taken to find impacts: " << elapsed.count() << "ms" << std::endl;
      std::cout << "Number of impacts: " << impactPoints.size() << std::endl;

      Eigen::MatrixXd P(impactPoints.size(), 3);
      result = Ot;
      for (int i = 0; i < impactPoints.size(); i++)
      {
        int idx = impactPoints[i].first;
        P.row(i) = Vs1.row(idx);

        if (m.corr(idx) != -1)
        {
          Eigen::Vector3d corrPoint = Ot.row(m.corr(idx));
          Eigen::Vector3d forceVec = (impactPoints[i].second * bs2 / brush_strength);

          applyKVL(
              result,
              corrPoint,
              forceVec,
              igl::KelvinletParams<double>(br2, scale, brushType),
              result);
        }
      }

      viewer.data_list[0].add_points(P, Eigen::RowVector3d(1, 0, 0));
      viewer.data_list[2].set_vertices(result);
      viewer.data_list[2].compute_normals();

      break;
    }
    default:
      std::cout << "Exiting menu" << std::endl;
      break;
    }
  };

  viewer.callback_init = [&](igl::opengl::glfw::Viewer &)
  {
    viewer.core().viewport = Eigen::Vector4f(0, 0, 300, 600); // Core 1 (left)
    core1 = viewer.core_list[0].id;
    core2 = viewer.append_core(Eigen::Vector4f(300, 0, 300, 600)); // Core 2 (center)
    core3 = viewer.append_core(Eigen::Vector4f(600, 0, 300, 600)); // Core 3 (right)

    viewer.data_list[1].set_visible(false, core1);
    viewer.data_list[2].set_visible(false, core1);
    viewer.data_list[3].set_visible(false, core1);
    viewer.data_list[0].set_visible(false, core2);
    viewer.data_list[2].set_visible(false, core2);
    viewer.data_list[3].set_visible(false, core2);
    viewer.data_list[0].set_visible(false, core3);
    viewer.data_list[1].set_visible(false, core3);
    viewer.data_list[3].set_visible(false, core3);
    return false;
  };

  viewer.callback_key_down = [&](igl::opengl::glfw::Viewer &, unsigned int key, int mod)
  {
    if (key == GLFW_KEY_M)
    {
      launchMenu();
      return true;
    }
    else if (key == GLFW_KEY_SPACE)
    {
      std::cout << "Switching mode: ";
      if (mode == MRK)
      {
        viewer.data_list[0].clear_points();
        viewer.data_list[2].clear_points();
        mode = CRS;
        std::cout << "MRK -> CRS" << std::endl;
      }
      else if (mode == CRS)
      {
        viewer.data_list[0].clear_points();
        viewer.data_list[2].clear_points();
        mode = DFM;
        std::cout << "CRS -> DFM" << std::endl;
      }
      else if (mode == DFM)
      {
        Vt = Ot;
        viewer.data_list[2].set_vertices(Vt);
        viewer.data_list[2].compute_normals();
        viewer.data_list[0].clear_points();
        viewer.data_list[2].clear_points();
        viewer.data_list[0].add_points(m.getFMarked(), Eigen::RowVector3d(0, 1, 0));
        viewer.data_list[2].add_points(m.getRMarked(), Eigen::RowVector3d(0, 1, 0));
        mode = MRK;
        std::cout << "DFM -> MRK" << std::endl;
      }
    }
    return false;
  };

  viewer.callback_post_resize = [&](igl::opengl::glfw::Viewer &v, int w, int h)
  {
    v.core(core1).viewport = Eigen::Vector4f(0, 0, w / 3, h);         // Left
    v.core(core2).viewport = Eigen::Vector4f(w / 3, 0, w / 3, h);     // Center
    v.core(core3).viewport = Eigen::Vector4f(2 * w / 3, 0, w / 3, h); // Right
    return true;
  };

  int v1 = -1, v2 = -1;
  viewer.callback_mouse_down =
      [&](igl::opengl::glfw::Viewer &viewer, int, int) -> bool
  {
    Eigen::Vector3f bc;
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;

    int fid;
    auto x = viewer.current_mouse_x;
    auto y =
        viewer.core().viewport(3) - static_cast<float>(viewer.current_mouse_y);

    if (viewer.core().id == core1)
    {
      V = Vs1;
      F = Fs;
    }
    else if (viewer.core().id == core2)
    {
      V = Vs2;
      F = Fs;
    }
    else if (viewer.core().id == core3)
    {
      V = Vt;
      F = Ft;
    }

    if (igl::unproject_onto_mesh(Eigen::Vector2f(x, y),
                                 viewer.core().view,
                                 viewer.core().proj,
                                 viewer.core().viewport,
                                 V,
                                 F,
                                 fid,
                                 bc))
    {
      if (viewer.core().id == core1 && mode == MRK)
      {
        int v = m.findClosestVertex(0, fid, bc);
        if (m.get(v) == -1)
        {
          if (v1 == v || v1 != -1)
          {
            viewer.data_list[0].clear_points();
            viewer.data_list[0].add_points(m.getFMarked(), Eigen::RowVector3d(0, 1, 0));
          }
          if (v1 != v)
          {
            Eigen::MatrixXd P(1, 3);
            P.row(0) = V.row(v);
            viewer.data_list[0].add_points(P, Eigen::RowVector3d(1, 0, 0));
            v1 = v;
            v2 = -1;
          }
          else
          {
            v1 = -1;
          }
        }
        else
        {
          m.erase(v);
          viewer.data_list[0].clear_points();
          viewer.data_list[0].add_points(m.getFMarked(), Eigen::RowVector3d(0, 1, 0));
          viewer.data_list[2].clear_points();
          viewer.data_list[2].add_points(m.getRMarked(), Eigen::RowVector3d(0, 1, 0));

          if (v1 == v)
          {
            v1 = -1;
            v2 = -1;
          }
        }

        return true;
      }
      else if (viewer.core().id == core3 && v1 != -1)
      {
        int v = m.findClosestVertex(1, fid, bc);
        m.set(v1, v);
        viewer.data_list[2].clear_points();
        viewer.data_list[2].add_points(m.getRMarked(), Eigen::RowVector3d(0, 1, 0));
        viewer.data_list[0].clear_points();
        viewer.data_list[0].add_points(m.getFMarked(), Eigen::RowVector3d(0, 1, 0)); // changing red to green

        return true;
      }
      else if (viewer.core().id == core1 && mode == CRS)
      {
        int v = m.findClosestVertex(0, fid, bc);
        int c = m.corr(v);

        if (c != -1)
        {
          viewer.data_list[0].clear_points();
          Eigen::MatrixXd P(1, 3);
          P.row(0) = V.row(v);
          viewer.data_list[0].add_points(P, Eigen::RowVector3d(0, 0, 1));

          viewer.data_list[2].clear_points();
          P.row(0) = Vt.row(c);
          viewer.data_list[2].add_points(P, Eigen::RowVector3d(0, 0, 1));
        }

        return true;
      }
      else if (viewer.core().id == core2 && mode == DFM)
      {
        projectedPoint = igl::unproject(
                             Eigen::Vector3f(x, y, viewer.down_mouse_z),
                             viewer.core().view,
                             viewer.core().proj,
                             viewer.core().viewport)
                             .template cast<double>();
        posStart = bc(0) * V.row(F(fid, 0)) +
                   bc(1) * V.row(F(fid, 1)) +
                   bc(2) * V.row(F(fid, 2));
      }
    }

    return false;
  };

  viewer.callback_mouse_move =
      [&](igl::opengl::glfw::Viewer &viewer, int, int) -> bool
  {
    if (!posStart.isZero() && !posStart.hasNaN())
    {
      posEnd = igl::unproject(
                   Eigen::Vector3f(viewer.current_mouse_x,
                                   viewer.core().viewport[3] -
                                       static_cast<float>(viewer.current_mouse_y),
                                   viewer.down_mouse_z),
                   viewer.core().view,
                   viewer.core().proj,
                   viewer.core().viewport)
                   .template cast<double>();
      posEnd = posStart + (posEnd - projectedPoint);

      // exaggerate the force by a little bit
      Eigen::Vector3d forceVec = (posEnd - posStart) * brush_strength;

      applyKVL(
          Vs2,
          posStart,
          forceVec,
          igl::KelvinletParams<double>(brushRadius, scale, brushType),
          result);
      viewer.data_list[1].set_vertices(result);
      viewer.data_list[1].compute_normals();
      return true;
    }
    return false;
  };

  viewer.callback_mouse_up =
      [&](igl::opengl::glfw::Viewer &viewer, int, int) -> bool
  {
    if (!posStart.isZero())
    {
      Vs2 = viewer.data_list[1].V;
      posStart.setZero();
      return true;
    }
    return false;
  };

  std::cout << "Press 'M' to open the menu" << std::endl;
  std::cout << "Press 'SPACE' to switch modes" << std::endl;
  std::cout << "Press 'ESC' to exit" << std::endl;

  viewer.launch(false, "Pose Transfer for 3D Objects", 900, 600);
  return EXIT_SUCCESS;
}

void applyKVL(
    Eigen::MatrixXd &V,
    Eigen::Vector3d &posStart,
    Eigen::Vector3d &forceVec,
    const igl::KelvinletParams<double> &params,
    Eigen::MatrixXd &result)
{
  Eigen::Matrix3d twist, pinch;
  twist << 0, 1, -1, -1, 0, 1, 1, -1, 0; // skew-symmetric
  pinch << 0, 1, 1, 1, 0, 1, 1, 1, 0;    // symmetric

  int scaleFactor = forceVec.norm();
  Eigen::Matrix3d mat;
  switch (params.brushType)
  {
  case igl::BrushType::GRAB:
    mat.setZero();
    break;
  case igl::BrushType::SCALE:
    mat = Eigen::Matrix3d::Identity() * scaleFactor;
    break;
  case igl::BrushType::TWIST:
    mat = twist * scaleFactor;
    break;
  case igl::BrushType::PINCH:
    mat = pinch * scaleFactor;
    break;
  }

  igl::kelvinlets(
      V,
      posStart,
      forceVec,
      mat,
      params,
      result);
}
