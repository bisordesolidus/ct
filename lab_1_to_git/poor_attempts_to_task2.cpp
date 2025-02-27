// -----------------------------------------------------------------------------
//
//  Gmsh C++ tutorial 13
//
//  Remeshing an STL file without an underlying CAD model
//
// -----------------------------------------------------------------------------

#include <set>
#include <cmath>
#include <gmsh.h>

int main(int argc, char **argv)
{
  gmsh::initialize();

  gmsh::model::add("t13");

  // Let's merge an STL mesh that we would like to remesh (from the parent
  // directory):
  try {
    gmsh::merge("../tarakan.STL");
  } catch(...) {
    gmsh::logger::write("Could not load STL mesh: bye!");
    gmsh::finalize();
    return 0;
  }

  // We first classify ("color") the surfaces by splitting the original surface
  // along sharp geometrical features. This will create new discrete surfaces,
  // curves and points.


  // Angle between two triangles above which an edge is considered as sharp:
  double angle = 40;

  // For complex geometries, patches can be too complex, too elongated or too
  // large to be parametrized; setting the following option will force the
  // creation of patches that are amenable to reparametrization:
  bool forceParametrizablePatches = true;

  // For open surfaces include the boundary edges in the classification process:
  bool includeBoundary = true;

  // Force curves to be split on given angle:
  double curveAngle = 180;

  gmsh::model::mesh::classifySurfaces(angle * M_PI / 180., includeBoundary,
                                      forceParametrizablePatches,
                                      curveAngle * M_PI / 180.);

  // Create a geometry for all the discrete curves and surfaces in the mesh, by
  // computing a parametrization for each one
  gmsh::model::mesh::createGeometry();

  // Create a volume from all the surfaces
  
  double mesh_size = 0.4;
  gmsh::option::setNumber("Mesh.MeshSizeMax", mesh_size);

  std::vector<std::pair<int, int> > s;
  gmsh::model::getEntities(s, 2);
  gmsh::model::mesh::setSize(s, mesh_size);
  std::vector<int> sl;
  for(auto surf : s) sl.push_back(surf.second);
  int l = gmsh::model::geo::addSurfaceLoop(sl);
  gmsh::model::geo::addVolume({l});

  gmsh::model::geo::synchronize();

  // We specify element sizes imposed by a size field, just because we can :-)
  bool funny = false; // false;
  int f = gmsh::model::mesh::field::add("MathEval");
  if(funny)
    gmsh::model::mesh::field::setString(f, "F", "2*Sin((x+y)/5) + 3");
  else
    gmsh::model::mesh::field::setString(f, "F", "4");
  gmsh::model::mesh::field::setAsBackgroundMesh(f);

  gmsh::model::mesh::generate(3);

  gmsh::write("hand.msh");

  // Launch the GUI to see the results:
  std::set<std::string> args(argv, argv + argc);
  if(!args.count("-nopopup")) gmsh::fltk::run();

  
  gmsh::finalize();
  return 0;
}




























/*#include <set>
#include <cmath>
#include <gmsh.h>
#include <iostream>


int main(int argc, char** argv) {
    gmsh::initialize();
    gmsh::option::setNumber("General.Terminal", 1);
    gmsh::model::add("iron_man");
    //gmsh::option::setNumber("Mesh.CharacteristicLengthMin", 0.1);
    //gmsh::option::setNumber("Mesh.CharacteristicLengthMax", 0.4);
    gmsh::option::setNumber("Mesh.MeshSizeMax", 1);
    std::string stl_file = "../merge_file.stl";
    try {
        gmsh::merge(stl_file);
    } catch (const std::exception& e) {
        std::cerr << "Error" << std::endl;
        gmsh::finalize();
        return 1;
    }



    //gmsh::option::setNumber("Mesh.AngleToleranceFacetOverlap", 5.0);
    //gmsh::option::setNumber("Mesh.Tolerance", 0.01);
    // Создаем топологию
  //gmsh::model::mesh::createTopology();
    //gmsh::model::mesh::createGeometry();
  double angle = 60;

  bool forceParametrizablePatches = true;

  bool includeBoundary = true;
  double curveAngle = 180;

  gmsh::model::mesh::classifySurfaces(angle * M_PI / 180., includeBoundary,
                                      forceParametrizablePatches,
                                      curveAngle * M_PI / 180.);
  gmsh::model::mesh::createGeometry();


  std::vector<std::pair<int, int> > s;
  gmsh::model::getEntities(s, 2);
  std::vector<int> sl;
  for(auto surf : s) sl.push_back(surf.second);
  int l = gmsh::model::geo::addSurfaceLoop(sl);
  gmsh::model::geo::addVolume({l});

  gmsh::model::geo::synchronize();
  gmsh::model::mesh::generate(3);

    gmsh::write("neko.msh");

    gmsh::finalize();
    return 0;
}*/
    /*gmsh::initialize();
    gmsh::option::setNumber("General.Terminal", 1);

    std::string stl_file = "../merge_file_fixed.stl"; 
    try {
        gmsh::merge(stl_file);
    } catch (const std::exception& e) {
        std::cerr << "Error merging STL file: " << e.what();
        gmsh::finalize();
        return 1;
    }

    gmsh::option::setNumber("Mesh.CharacteristicLengthMin", 0.05); // Минимальный размер элемента
    gmsh::option::setNumber("Mesh.CharacteristicLengthMax", 1.0); // Максимальный размер элемента
    bool includeBoundary = true;

    double curveAngle = 180;
    double angle = 40;
    bool forceParametrizablePatches = false;
    gmsh::model::mesh::classifySurfaces(angle * M_PI / 180., includeBoundary,
                                      forceParametrizablePatches,
                                      curveAngle * M_PI / 180.);
    
    gmsh::option::setNumber("Mesh.Algorithm", 5);
    gmsh::model::mesh::createGeometry();

    std::vector<std::pair<int, int> > s;
    gmsh::model::getEntities(s, 2);
    std::vector<int> sl;
    for(auto surf : s) sl.push_back(surf.second);
    int l = gmsh::model::geo::addSurfaceLoop(sl);
    int f = gmsh::model::mesh::field::add("MathEval");
    gmsh::model::geo::addVolume({l});
    gmsh::model::geo::synchronize();
    gmsh::model::mesh::field::setString(f, "F", "4");
    gmsh::model::mesh::field::setAsBackgroundMesh(f);
    gmsh::model::mesh::generate(3);
    gmsh::write("neko.msh");


    gmsh::finalize();
    return 0;*/


/*int main(int argc, char **argv)
{
  gmsh::initialize();

  gmsh::model::add("t13");

  // Let's merge an STL mesh that we would like to remesh (from the parent
  // directory):
  try {
    gmsh::merge("../merge_file.stl");
  } catch(...) {
    gmsh::logger::write("Could not load STL mesh: bye!");
    gmsh::finalize();
    return 0;
  }

  // We first classify ("color") the surfaces by splitting the original surface
  // along sharp geometrical features. This will create new discrete surfaces,
  // curves and points.


  // Angle between two triangles above which an edge is considered as sharp:
  double angle = 40;

  // For complex geometries, patches can be too complex, too elongated or too
  // large to be parametrized; setting the following option will force the
  // creation of patches that are amenable to reparametrization:
  bool forceParametrizablePatches = false;

  // For open surfaces include the boundary edges in the classification process:
  bool includeBoundary = true;

  // Force curves to be split on given angle:
  double curveAngle = 180;

  gmsh::model::mesh::classifySurfaces(angle * M_PI / 180., includeBoundary,
                                      forceParametrizablePatches,
                                      curveAngle * M_PI / 180.);

  // Create a geometry for all the discrete curves and surfaces in the mesh, by
  // computing a parametrization for each one
  gmsh::model::mesh::createGeometry();

  // Create a volume from all the surfaces
  std::vector<std::pair<int, int> > s;
  gmsh::model::getEntities(s, 2);
  std::vector<int> sl;
  for(auto surf : s) sl.push_back(surf.second);
  int l = gmsh::model::geo::addSurfaceLoop(sl);
  gmsh::model::geo::addVolume({l});

  gmsh::model::geo::synchronize();

  // We specify element sizes imposed by a size field, just because we can :-)
  bool funny = true; // false;
  int f = gmsh::model::mesh::field::add("MathEval");
  if(funny)
    gmsh::model::mesh::field::setString(f, "F", "2*Sin((x+y)/5) + 3");
  else
    gmsh::model::mesh::field::setString(f, "F", "4");
  gmsh::model::mesh::field::setAsBackgroundMesh(f);

  gmsh::model::mesh::generate(3);

  gmsh::write("t5.msh");

  // Launch the GUI to see the results:
  std::set<std::string> args(argv, argv + argc);
  if(!args.count("-nopopup")) gmsh::fltk::run();

  gmsh::finalize();
  return 0;
}*/

