#include <iostream>
#include <fstream>
#include <string>
#include <gmsh.h>
#include <vector>

int main(int argc, char **argv) {
    gmsh::initialize();
    double R = 10;
    double r = 7;
    double rr = 5; 
    int legs = gmsh::model::occ::addDisk(0,0,0,R ,R);
    int body = gmsh::model::occ::addDisk((R+2*r+rr)-2.6,0,0,rr ,rr);
    int head = gmsh::model::occ::addDisk(R+r-1.6,0,0,r ,r);
    gmsh::option::setNumber("Mesh.CharacteristicLengthMin", 0.05);
    gmsh::option::setNumber("Mesh.CharacteristicLengthMax", 0.1);
    gmsh::vectorpair objectDimTags = {{2, static_cast<int>(legs)}};
    gmsh::vectorpair toolDimTags = {{2, static_cast<int>(body)}, {2, static_cast<int>(head)}};
    gmsh::vectorpair outDimTags;
    std::vector<gmsh::vectorpair> outDimTagsMap;
    bool removeObject = true;  // Удалить circle1 после объединения
    bool removeTool   = true;  // Удалить circle2 и circle3 после объединения
    int tag = -1; //auto generate tag
    gmsh::model::occ::fuse(objectDimTags, toolDimTags, outDimTags, outDimTagsMap, tag, removeObject, removeTool);

    gmsh::model::occ::synchronize();
    gmsh::model::mesh::generate(2);
    gmsh::write("figure.msh");

    gmsh::finalize();
    return 0;



}  
    /*gmsh::initialize(argc, argv);
    gmsh::model::add("square_with_hole");

    double lc = 0.1;       // Characteristic length
    double squareSize = 1.0; // Size of the square
    double circleRadius = 0.3; // Radius of the circle
    double circleCenterX = squareSize / 2.0;
    double circleCenterY = squareSize / 2.0;

    // 1. Создание квадрата (geo)
    int p1 = gmsh::model::geo::addPoint(0, 0, 0, lc);
    int p2 = gmsh::model::geo::addPoint(squareSize, 0, 0, lc);
    int p3 = gmsh::model::geo::addPoint(squareSize, squareSize, 0, lc);
    int p4 = gmsh::model::geo::addPoint(0, squareSize, 0, lc);

    int l1 = gmsh::model::geo::addLine(p1, p2);
    int l2 = gmsh::model::geo::addLine(p2, p3);
    int l3 = gmsh::model::geo::addLine(p3, p4);
    int l4 = gmsh::model::geo::addLine(p4, p1);

    int curveLoopSquare = gmsh::model::geo::addCurveLoop({l1, l2, l3, l4});
    int squareSurface = gmsh::model::geo::addPlaneSurface({curveLoopSquare});

    // 2. Создание окружности (occ)
    int circle = gmsh::model::occ::addDisk(circleCenterX, circleCenterY, 0, circleRadius, circleRadius);

    // 3. Fragment
    std::vector<std::pair<int, int>> objectToFragment = {{2, static_cast<int>(squareSurface)}}; // Square (2D Entity)
    std::vector<std::pair<int, int>> toolToFragment = {{2, static_cast<int>(circle)}}; // Circle (OCC Entity)

    std::vector<std::vector<std::size_t>> outMasterTags;
    std::vector<std::pair<int, int>> out;
    std::vector<std::vector<std::size_t>> outtoolTags;
    std::vector<std::vector<std::pair<int, int>>> outDimTags; 
    int removeObject = -1;
    bool removeTool = true;  //true=удалять, false=не удалять
    bool frackSurface = false; // Разбивать ли поверхности

    gmsh::model::occ::fragment(objectToFragment, toolToFragment, out, outDimTags, removeObject, removeTool, frackSurface);

    // 4. Synchronize
    gmsh::model::occ::synchronize();

    //5. Add Physical Groups (VERY IMPORTANT!)
    // Find all surfaces after fragmentation
    std::vector<std::pair<int, int>> allSurfaces;
    gmsh::model::getEntities(allSurfaces, 2);

    std::vector<int> surfaceTags;
    for (const auto& tag : allSurfaces) {
        surfaceTags.push_back(tag.second);
    }

    // Create a physical surface for the domain (the square with the hole)
    gmsh::model::addPhysicalGroup(2, surfaceTags, 1); //Surface
    gmsh::model::setPhysicalName(2, 1, "Domain");

    // Create physical lines for boundaries
    std::vector<std::pair<int, int>> allCurves;
    gmsh::model::getEntities(allCurves, 1);

    std::vector<int> curveTags;
    for (const auto& tag : allCurves) {
        curveTags.push_back(tag.second);
    }
    gmsh::model::addPhysicalGroup(1, curveTags, 2); //Curve
    gmsh::model::setPhysicalName(1, 2, "Boundary");

    // 6. Mesh generation
    gmsh::model::mesh::generate(2);

    // 7. Write to file
    gmsh::write("square_with_hole.msh");

    // 8. Visualization (optional)
    if (argc > 1 && std::string(argv[1]) == "-gui") {
        gmsh::fltk::run();
    }

    gmsh::finalize();
    return 0;
}*/


    /*double a = 5; //square size
    double r = 2;
    double lc = 0.01;
    gmsh::initialize();
    gmsh::model::add("figure");
    int p1 = gmsh::model::geo::addPoint(0, 0, 0, lc);
    int p2 = gmsh::model::geo::addPoint(a, 0, 0, lc);
    int p3 = gmsh::model::geo::addPoint(a, a, 0, lc);
    int p4 = gmsh::model::geo::addPoint(0, a, 0, lc);

    int l1 = gmsh::model::geo::addLine(p1, p2);
    int l2 = gmsh::model::geo::addLine(p2, p3);
    int l3 = gmsh::model::geo::addLine(p3, p4);
    int l4 = gmsh::model::geo::addLine(p4, p1);

    int curveLoop = gmsh::model::geo::addCurveLoop({l1, l2, l3, l4});
    int tobecut = gmsh::model::geo::addPlaneSurface({curveLoop});
    gmsh::model::geo::synchronize();

    int cutting = gmsh::model::occ::addDisk(r,r,0,r ,r);
    std::vector<std::pair<int, int>> objects = {{2, tobecut}};
    std::vector<std::pair<int, int>> tools = {{2, cutting}};
    std::vector<std::pair<int, int>> out;
    std::vector<std::vector<std::pair<int, int>>> outDimTags; 
    int removeObject = -1;
    bool removeTool = true;  //true=удалять, false=не удалять
    bool frackSurface = false; // Разбивать ли поверхности

    gmsh::model::occ::fragment(objects, tools, out, outDimTags, removeObject, removeTool, frackSurface);
    gmsh::model::occ::synchronize();
    gmsh::model::mesh::generate(2);
    gmsh::write("figure.msh");

    gmsh::finalize();
    return 0;*/