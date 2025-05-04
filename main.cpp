#define IGL_VIEWER_VIEWER_QUIET 0
#include <igl/opengl/glfw/Viewer.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/unproject.h>
#include <igl/kelvinlets.h>
#include <iostream>
#include <chrono>
#include "dfm_breakdown.h"

unsigned left_view, right_view;
typedef enum
{
  DFM,
  TRF
} OperatingMode;

void applyKVL(
    Eigen::MatrixXd &V,
    Eigen::Vector3d &posStart,
    Eigen::Vector3d &forceVec,
    const igl::KelvinletParams<double> &params,
    Eigen::MatrixXd &result);

int main(int argc, char *argv[])
{
  if (argc < 2)
  {
    std::cout << "Specify valid mesh file" << std::endl;
    exit(1);
  }

  igl::opengl::glfw::Viewer viewer;
  Eigen::MatrixXd Ov, V1, V2, V3;
  Eigen::MatrixXi F1, F2, F3;

  auto brushRadius = 1.;
  auto brushType = igl::BrushType::GRAB;
  auto brush_strength = 1.;
  auto scale = 1;
  OperatingMode mode = DFM;

  Eigen::Vector3d posStart(0, 0, 0), posEnd, uPoint, ps2;
  decltype(V1) result;

  viewer.load_mesh_from_file(argv[1]);
  viewer.load_mesh_from_file(argv[1]);
  // viewer.load_mesh_from_file(argv[1]);

  V1 = viewer.data_list[0].V;
  F1 = viewer.data_list[0].F;
  V2 = viewer.data_list[1].V;
  F2 = viewer.data_list[1].F;
  // V3 =
  Ov = V1;

  brush_strength = (V1.colwise().maxCoeff() - V1.colwise().minCoeff()).norm();
  brushRadius = brush_strength / 20;

  const auto launchMenu = [&]()
  {
    std::cout << "Menu" << std::endl;
    std::cout << "1. Change brush radius" << std::endl;
    std::cout << "2. Change brush strength" << std::endl;
    std::cout << "3. Change scale" << std::endl;
    std::cout << "4. Reset deformation" << std::endl;
    std::cout << "5. Transfer deformation" << std::endl;
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
      V1 = Ov;
      V2 = Ov;
      viewer.data_list[0].set_vertices(V1);
      viewer.data_list[1].set_vertices(V2);
      viewer.data_list[0].compute_normals();
      viewer.data_list[1].compute_normals();
      break;
    case 5:
    {
      // mode = TRF;
      // std::cout << "Click on point in right mesh" << std::endl;

      // Eigen::MatrixXd P = impactPoints(Ov, V1);
      // std::cout << "Number of vertices changed: " << P.rows() << std::endl;
      // viewer.data_list[1].add_points(P, Eigen::RowVector3d(1, 0, 0));

      auto start = std::chrono::high_resolution_clock::now();
      auto out = find_impacts_and_forces(
          Ov,
          V1, 0.25);
      auto end = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double, std::milli> elapsed = end - start;

      std::cout << "Impacts on " << out.size() << " vertices" << std::endl;
      result = Ov;
      Eigen::MatrixXd t;
      for (int i = 0; i < out.size(); i++)
      {
        Eigen::Vector3d point = out[i].first.transpose();
        applyKVL(
            result,
            point,
            out[i].second,
            igl::KelvinletParams<double>(brushRadius, scale, brushType),
            t);
        result = t;
      }

      viewer.data_list[1].set_vertices(result);
      viewer.data_list[1].compute_normals();
      std::cout << "DT complete" << std::endl;
      std::cout << elapsed.count() << ",";
      printError(V1, result);

      break;
    }
    default:
      std::cout << "Exiting menu" << std::endl;
      break;
    }
  };

  viewer.callback_init = [&](igl::opengl::glfw::Viewer &)
  {
    viewer.core().viewport = Eigen::Vector4f(0, 0, 640, 800);
    left_view = viewer.core_list[0].id;
    right_view = viewer.append_core(Eigen::Vector4f(640, 0, 640, 800));
    // By default, when a core is appended, all loaded meshes will be displayed in that core.
    viewer.data_list[0].set_visible(false, right_view);
    viewer.data_list[1].set_visible(false, left_view);
    return false;
  };

  viewer.callback_key_down = [&](igl::opengl::glfw::Viewer &, unsigned int key, int mod)
  {
    if (key == GLFW_KEY_M)
    {
      launchMenu();
      return true;
    }
    return false;
  };

  viewer.callback_post_resize = [&](igl::opengl::glfw::Viewer &v, int w, int h)
  {
    v.core(left_view).viewport = Eigen::Vector4f(0, 0, w / 2, h);
    v.core(right_view).viewport = Eigen::Vector4f(w / 2, 0, w - (w / 2), h);
    return true;
  };

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

    if (viewer.core().id == left_view)
    {
      V = V1;
      F = F1;
    }
    else
    {
      V = V2;
      F = F2;
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

      if (viewer.core().id == left_view && mode == DFM)
      {
        uPoint = igl::unproject(Eigen::Vector3f(x, y, viewer.down_mouse_z),
                                viewer.core().view,
                                viewer.core().proj,
                                viewer.core().viewport)
                     .template cast<double>();
        Eigen::Vector3d clicked_point = bc(0) * V.row(F(fid, 0)) +
                                        bc(1) * V.row(F(fid, 1)) +
                                        bc(2) * V.row(F(fid, 2));
        posStart = clicked_point;

        return true;
      }
      else if (viewer.core().id == right_view && mode == TRF)
      {
        uPoint = igl::unproject(Eigen::Vector3f(x, y, viewer.down_mouse_z),
                                viewer.core().view,
                                viewer.core().proj,
                                viewer.core().viewport)
                     .template cast<double>();
        Eigen::Vector3d clicked_point = bc(0) * V.row(F(fid, 0)) +
                                        bc(1) * V.row(F(fid, 1)) +
                                        bc(2) * V.row(F(fid, 2));

        auto start = std::chrono::high_resolution_clock::now();
        Eigen::Vector3d force = estimate_force_vector(Ov, V1, clicked_point.transpose(), 0.01);
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> elapsed = end - start;
        std::cout << "Recovered force: " << force.transpose() << "\nTime taken: " << elapsed.count() << std::endl;

        applyKVL(
            V2,
            clicked_point,
            force,
            igl::KelvinletParams<double>(brushRadius, scale, brushType),
            result);
        viewer.data_list[1].set_vertices(result);
        viewer.data_list[1].compute_normals();
        mode = DFM;
      }
    }

    return false;
  };

  viewer.callback_mouse_move =
      [&](igl::opengl::glfw::Viewer &viewer, int, int) -> bool
  {
    if (!posStart.isZero() && !posStart.hasNaN() && viewer.core().id == left_view)
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
      posEnd = posStart + (posEnd - uPoint);

      // exaggerate the force by a little bit
      Eigen::Vector3d forceVec = (posEnd - posStart) * brush_strength;

      applyKVL(
          V1,
          posStart,
          forceVec,
          igl::KelvinletParams<double>(brushRadius, scale, brushType),
          result);
      viewer.data_list[0].set_vertices(result);

      viewer.data_list[0].compute_normals();
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
      V1 = viewer.data_list[0].V;
      posStart.setZero();
      return true;
    }
    return false;
  };

  viewer.launch();
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
