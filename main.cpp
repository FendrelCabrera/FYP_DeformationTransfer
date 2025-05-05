#define IGL_VIEWER_VIEWER_QUIET 0
#include <igl/opengl/glfw/Viewer.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/unproject.h>
#include <igl/kelvinlets.h>
#include <iostream>
#include <chrono>
#include "dfm_breakdown.h"
#include "mapper.cpp"

unsigned core1, core2, core3, core4;

void applyKVL(
    Eigen::MatrixXd &V,
    Eigen::Vector3d &posStart,
    Eigen::Vector3d &forceVec,
    const igl::KelvinletParams<double> &params,
    Eigen::MatrixXd &result);

int main(int argc, char *argv[])
{
  if (argc < 4)
  {
    std::cout << "Specify valid mesh files" << std::endl;
    exit(1);
  }

  igl::opengl::glfw::Viewer viewer;
  Eigen::MatrixXd V1, V2, V3, V4;
  Eigen::MatrixXi F12, F34;

  auto brushRadius = 1.;
  auto brushType = igl::BrushType::GRAB;
  auto brush_strength = 1.;
  auto scale = 1;

  Eigen::Vector3d posStart(0, 0, 0), posEnd, uPoint, ps2;
  decltype(V1) result;

  viewer.load_mesh_from_file(argv[1]);
  viewer.load_mesh_from_file(argv[2]);
  viewer.load_mesh_from_file(argv[3]);
  viewer.load_mesh_from_file(argv[3]);

  V1 = viewer.data_list[0].V;
  F12 = viewer.data_list[0].F;
  V2 = viewer.data_list[1].V;
  V3 = V4 = viewer.data_list[2].V;
  F34 = viewer.data_list[2].F;

  std::cout << "(" << V1.rows() << "," << F12.rows() << ") (" << V3.rows() << "," << F34.rows() << ")" << std::endl;

  Mapper m(V1, F12, V3, F34);
  auto start = std::chrono::high_resolution_clock::now();
  auto impactPoints = find_impacts_and_forces(
      V1,
      V2, 0.25);
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double, std::milli> elapsed = end - start;
  std::cout << "Time taken to find impacts: " << elapsed.count() << "ms" << std::endl;
  std::cout << "Number of impacts: " << impactPoints.size() << std::endl;

  brush_strength = (V1.colwise().maxCoeff() - V1.colwise().minCoeff()).norm();
  brushRadius = brush_strength / 20;

  const auto launchMenu = [&]()
  {
    std::cout << "Menu" << std::endl;
    std::cout << "1. Change brush radius" << std::endl;
    std::cout << "2. Change brush strength" << std::endl;
    std::cout << "3. Change scale" << std::endl;
    std::cout << "4. Reset deformation" << std::endl;
    std::cout << "5. Apply deformation" << std::endl;
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
      V4 = V3;
      viewer.data_list[3].set_vertices(V4);
      viewer.data_list[3].compute_normals();
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

      break;
    }
    default:
      std::cout << "Exiting menu" << std::endl;
      break;
    }
  };

  viewer.callback_init = [&](igl::opengl::glfw::Viewer &)
  {
    viewer.core().viewport = Eigen::Vector4f(0, 400, 640, 400); // Core 1 (top-left)
    core1 = viewer.core_list[0].id;
    core2 = viewer.append_core(Eigen::Vector4f(640, 400, 640, 400)); // Core 2 (top-right)
    core3 = viewer.append_core(Eigen::Vector4f(0, 0, 640, 400));     // Core 3 (bottom-left)
    core4 = viewer.append_core(Eigen::Vector4f(640, 0, 640, 400));   // Core 4 (bottom-right)

    viewer.data_list[0].set_visible(false, core2);
    viewer.data_list[0].set_visible(false, core3);
    viewer.data_list[0].set_visible(false, core4);
    viewer.data_list[1].set_visible(false, core1);
    viewer.data_list[1].set_visible(false, core3);
    viewer.data_list[1].set_visible(false, core4);
    viewer.data_list[2].set_visible(false, core1);
    viewer.data_list[2].set_visible(false, core2);
    viewer.data_list[2].set_visible(false, core4);
    viewer.data_list[3].set_visible(false, core1);
    viewer.data_list[3].set_visible(false, core2);
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
    return false;
  };

  viewer.callback_post_resize = [&](igl::opengl::glfw::Viewer &v, int w, int h)
  {
    v.core(core1).viewport = Eigen::Vector4f(0, h / 2, w / 2, h - (h / 2));           // Top-left
    v.core(core2).viewport = Eigen::Vector4f(w / 2, h / 2, w - (w / 2), h - (h / 2)); // Top-right
    v.core(core3).viewport = Eigen::Vector4f(0, 0, w / 2, h / 2);                     // Bottom-left
    v.core(core4).viewport = Eigen::Vector4f(w / 2, 0, w - (w / 2), h / 2);           // Bottom-right
    return true;
  };

  int v1 = -1, v2 = -1;
  viewer.callback_mouse_down =
      [&](igl::opengl::glfw::Viewer &viewer, int, int) -> bool
  {
    std::cout << "Mouse down" << std::endl;
    Eigen::Vector3f bc;
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;

    int fid;
    auto x = viewer.current_mouse_x;
    auto y =
        viewer.core().viewport(3) - static_cast<float>(viewer.current_mouse_y);

    if (viewer.core().id == core1)
    {
      std::cout << "Click on point in left mesh" << std::endl;
      V = V1;
      F = F12;
    }
    else if (viewer.core().id == core2)
    {
      return false;
    }
    else if (viewer.core().id == core3)
    {
      std::cout << "Click on point in right mesh" << std::endl;
      V = V3;
      F = F34;
    }
    else if (viewer.core().id == core4)
    {
      return false;
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
      std::cout << "Projected to mesh" << std::endl;
      if (viewer.core().id == core1)
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
            std::cout << "Adding point" << P << std::endl;
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
        }

        return true;
      }
      else if (viewer.core().id == core3 && v1 != -1)
      {
        int v = m.findClosestVertex(1, fid, bc);
        if (v2 != -1 && v2 != v)
        {
          viewer.data_list[1].clear_points();
          viewer.data_list[1].add_points(m.getRMarked(), Eigen::RowVector3d(0, 1, 0));
        }
        if (v2 != v)
        {
          Eigen::MatrixXd P(1, 3);
          P.row(0) = V.row(v);
          viewer.data_list[1].add_points(P, Eigen::RowVector3d(0, 1, 0));
          v2 = v;
        }
        m.set(v1, v2);

        return true;
      }
    }

    return false;
  };

  viewer.callback_mouse_move =
      [&](igl::opengl::glfw::Viewer &viewer, int, int) -> bool
  {
    // if (!posStart.isZero() && !posStart.hasNaN() && viewer.core().id == left_view)
    // {
    //   posEnd = igl::unproject(
    //                Eigen::Vector3f(viewer.current_mouse_x,
    //                                viewer.core().viewport[3] -
    //                                    static_cast<float>(viewer.current_mouse_y),
    //                                viewer.down_mouse_z),
    //                viewer.core().view,
    //                viewer.core().proj,
    //                viewer.core().viewport)
    //                .template cast<double>();
    //   posEnd = posStart + (posEnd - uPoint);

    //   // exaggerate the force by a little bit
    //   Eigen::Vector3d forceVec = (posEnd - posStart) * brush_strength;

    //   applyKVL(
    //       V1,
    //       posStart,
    //       forceVec,
    //       igl::KelvinletParams<double>(brushRadius, scale, brushType),
    //       result);
    //   viewer.data_list[0].set_vertices(result);

    //   viewer.data_list[0].compute_normals();
    //   viewer.data_list[1].compute_normals();
    //   return true;
    // }
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
