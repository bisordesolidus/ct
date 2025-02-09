#include <iostream>
#include <gmsh.h>
#include <cmath>

/*int draw_circle(int r, int i, int x0 = 10, int y0 = 0, int z0 = 0,int lc = 0.1){
    gmsh::model::geo::addPoint(x0 + r, y0, z0, lc, (i-1)*4+1); // x = r, z = 0
    gmsh::model::geo::addPoint(x0, y0, z0+r, lc, (i-1)*4+2); // x = 0, z = r
    gmsh::model::geo::addPoint(x0 - r, y0, z0, lc, (i-1)*4+3); // x = -r, z = 0
    gmsh::model::geo::addPoint(x0, y0, z0-r, lc, (i-1)*4+4); // x = 0, z = -r

    gmsh::model::geo::addCircleArc(1, 10, 2, (i-1)*4+1); // Дуга от точки 1 до точки 2 через центр (0,0,0)
    gmsh::model::geo::addCircleArc(2, 10, 3, (i-1)*4+2); // Дуга от точки 2 до точки 3 через центр (0,0,0)
    gmsh::model::geo::addCircleArc(3, 10, 4, (i-1)*4+3); // Дуга от точки 3 до точки 4 через центр (0,0,0)
    gmsh::model::geo::addCircleArc(4, 10, 1, (i-1)*4+4); // Дуга от точки 4 до точки 1 через центр (0,0,0)

    std::vector<int> curveTags = {(i-1)*4+1, (i-1)*4+2, (i-1)*4+3, (i-1)*4+4};
    gmsh::model::geo::addCurveLoop(curveTags, i);
    return i;
}*/

int main(int argc, char **argv) {
    double lc = 0.1;           // Характерный размер элемента сетки
    double R = 10;            // Большой радиус тора
    double r2 = 7;
    double r1 = 3;            // Малый радиус трубки тора

    double x0 = 10;
    double y0 = 0;
    double z0 = 0;
    gmsh::initialize();
    gmsh::model::geo::addPoint(x0, y0, z0, lc, 10); //center
    int torus2 = 3;
    int torus1 = 4;
    gmsh::model::occ::addTorus(0, 0, 0, R, r2,torus2);
    gmsh::model::occ::addTorus(0, 0, 0, R, r1,torus1);

    //std::vector<int> surfaceTags = {torus2,-torus1};
    //gmsh::model::occ::addPlaneSurface(surfaceTags, 5);

    std::vector<std::pair<int, int>> objects = {{3, torus2}};
    std::vector<std::pair<int, int>> tools = {{3, torus1}};
    std::vector<std::pair<int, int>> out;
    std::vector<std::vector<std::pair<int, int>>> outDimTags; 
    int removeObject = -1;
    bool removeTool = true;  //true=удалять, false=не удалять
    bool frackSurface = false; // Разбивать ли поверхности

    gmsh::model::occ::fragment(objects, tools, out, outDimTags, removeObject, removeTool, frackSurface);


    /*std::vector<std::pair<int,int>> tobecut = {{3, torus2}};
    std::vector<std::pair<int,int>> tocut = {{3, torus1}};
    std::vector<std::pair<int,int>> result;
    std::vector<std::pair<int>> empty_tag_vector;
    std::vector<std::pair<double>> empty_angle_vector;
    gmsh::model::occ::cut(tobecut, tocut, result, empty_tag_vector, empty_angle_vector, true, true);*/

    gmsh::model::occ::synchronize();
    gmsh::model::mesh::generate(3);
    gmsh::write("torus_occ.msh");

    gmsh::finalize();
    return 0;

}