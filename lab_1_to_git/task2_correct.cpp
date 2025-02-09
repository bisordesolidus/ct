#include <iostream>
#include <gmsh.h>

int main(int argc, char **argv) {
    gmsh::initialize(argc, argv);
    gmsh::option::setNumber("General.Terminal", 1);

    // 1. Импорт STL файла
    std::string stl_file = "merge_file_fixed.stl";
    try {
        gmsh::merge(stl_file);
    } catch (const std::exception& e) {
        std::cerr << "Error merging STL file: " << e.what() << std::endl;
        gmsh::finalize();
        return 1;
    }

    gmsh::model::mesh::generate(3); 
    gmsh::write("niko.msh");
    gmsh::finalize();

    return 0;
}