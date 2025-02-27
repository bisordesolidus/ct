#include <iostream>
#include <cmath>
#include <vector>
#include <set>
#include <algorithm>

#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkTetra.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkSmartPointer.h>

#include <gmsh.h>

using namespace std;

namespace m{
        // ----------------------- Параметры модели -----------------------
        double scale = 1.0; // Масштаб всей модели
        //double head_height =  * scale;
        const double head_size = 6 * scale;
        const double body_width = 4 * scale;
        const double body_height = 7 * scale;
        const double arm_width = 3 * scale;
        const double arm_height = body_height;
        const double leg_width = 2 * scale;
        const double leg_height = 8 * scale;
        const double z_offset = 0;
    
}


bool is_head(std::vector<double> nodesCoord, unsigned i){
    double x = nodesCoord[3*(i-1)];
    double y = nodesCoord[3*(i-1)+1];
    double z = nodesCoord[3*(i-1)+2];
    if(x>=0 and x<=6 and z>=0 and z<=6 and y>0 and y<6){return true;}
    else{return false;}
}
bool is_body(std::vector<double> nodesCoord, unsigned i){
    double x = nodesCoord[3*(i-1)];
    double y = nodesCoord[3*(i-1)+1];
    double z = nodesCoord[3*(i-1)+2];
    if(x>=1 and x<=5 and z>=1 and z<=5 and y>6 and y<13){return true;}
    else{return false;}
}
bool is_right_hand(std::vector<double> nodesCoord, unsigned i){
    double x = nodesCoord[3*(i-1)];
    double y = nodesCoord[3*(i-1)+1];
    double z = nodesCoord[3*(i-1)+2];
    if(x>=6 and x<=8 and z>=1.5 and z<=4.5 and y>=6 and y<=13){return true;}
    else{return false;}
}
bool is_left_hand(std::vector<double> nodesCoord, unsigned i){
    double x = nodesCoord[3*(i-1)];
    double y = nodesCoord[3*(i-1)+1];
    double z = nodesCoord[3*(i-1)+2];
    if(x>=-2 and x<=1 and z>=1.5 and z<=4.5 and y>=6 and y<=13){return true;}
    else{return false;}
}
bool is_right_leg(std::vector<double> nodesCoord, unsigned i){
    double x = nodesCoord[3*(i-1)];
    double y = nodesCoord[3*(i-1)+1];
    double z = nodesCoord[3*(i-1)+2];
    if(x>=3 and x<=5 and z>=2 and z<=4 and y>=13 and y<=21){return true;}
    else{return false;}
}

bool is_left_leg(std::vector<double> nodesCoord, unsigned i){
    double x = nodesCoord[3*(i-1)];
    double y = nodesCoord[3*(i-1)+1];
    double z = nodesCoord[3*(i-1)+2];
    if(x>=1 and x<=3 and z>=2 and z<=4 and y>=13 and y<=21){return true;}
    else{return false;}
}


// Класс расчётной точки
class CalcNode
{
// Класс сетки будет friend-ом точки
friend class CalcMesh;

protected:
    // Координаты
    double x;
    double y;
    double z;
    // Некая величина, в попугаях
    double smth;
    // Скорость
    double vx;
    double vy;
    double vz;

    double wx;
    double wy;
    double wz;

public:
    // Конструктор по умолчанию
    CalcNode() : x(0.0), y(0.0), z(0.0), smth(0.0), vx(0.0), vy(0.0), vz(20.0)
    {
    }

    // Конструктор с указанием всех параметров
    CalcNode(double x, double y, double z, double smth, double vx, double vy, double vz =20) 
            : x(x), y(y), z(z), smth(smth), vx(vx), vy(vy), vz(vz)
    {
    }


    // Метод отвечает за перемещение точки
    // Движемся время tau из текущего положения с текущей скоростью
    void move(double tau, double parameter) {
        z += parameter * tau;
    }
    void rotate(double tau, int step, float wx, float ry, float rz){
        y+= -(wx*rz)*tau;
        z+=wx*ry*tau;
    }
};

// Класс элемента сетки
class Element
{
// Класс сетки будет friend-ом и элемента тоже
// (и вообще будет нагло считать его просто структурой)
friend class CalcMesh;

protected:
    // Индексы узлов, образующих этот элемент сетки
    unsigned long nodesIds[4];

};

// Класс расчётной сетки
class CalcMesh
{
protected:
    // 3D-сетка из расчётных точек
    vector<CalcNode> nodes;
    vector<Element> elements;


public:
    // Конструктор сетки из заданного stl-файла
    CalcMesh(const std::vector<double>& nodesCoords, const std::vector<std::size_t>& tetrsPoints) {

        // Пройдём по узлам в модели gmsh
        nodes.resize(nodesCoords.size() / 3);
        for(unsigned int i = 0; i < nodesCoords.size() / 3; i++) {
            // Координаты заберём из gmsh
            double pointX = nodesCoords[i*3];
            double pointY = nodesCoords[i*3 + 1];
            double pointZ = nodesCoords[i*3 + 2];
            // Модельная скалярная величина распределена как-то вот так
            double smth = pow(pointX, 2) + pow(pointY, 2) + pow(pointZ, 2);
            nodes[i] = CalcNode(pointX, pointY, pointZ, smth, 0.0, 0.0, 0.0);
        }

        // Пройдём по элементам в модели gmsh
        elements.resize(tetrsPoints.size() / 4);
        for(unsigned int i = 0; i < tetrsPoints.size() / 4; i++) {
            elements[i].nodesIds[0] = tetrsPoints[i*4] - 1;
            elements[i].nodesIds[1] = tetrsPoints[i*4 + 1] - 1;
            elements[i].nodesIds[2] = tetrsPoints[i*4 + 2] - 1;
            elements[i].nodesIds[3] = tetrsPoints[i*4 + 3] - 1;
        }
    }

    // Метод отвечает за выполнение для всей сетки шага по времени величиной tau
    void doTimeStep(double tau, int step,  std::vector<std::size_t>& nodeTagsIn1, std::vector<std::size_t>& nodeTagsIn2, std::vector<std::size_t>& nodeTagsIn3, std::vector<std::size_t>& nodeTagsIn4, std::vector<std::size_t>& nodeTagsIn5) {
        //сначала двигаем правую отн х руку потом левую. каждая по 5 тиков
        //1 -голова
        //2-торс
        //3,4 -руки
        //5,6 -ноги
        for(unsigned int i = 0; i < nodes.size(); i++) {
            nodes[i].smth =10e21 * pow(sin(nodes[i].x * nodes[i].x + nodes[i].y * nodes[i].y + nodes[i].z * nodes[i].z), 2)* sin(10e21 * tau*step);
            nodes[i].move(tau, 20);
        }
        double w=2.5;

        if(step%20 <5 and step%20>0){
                for (unsigned int i = 0; i < nodeTagsIn2.size(); i++)   {nodes[nodeTagsIn2[i]].rotate(tau,step, w , m::head_size, m::z_offset+ (m::head_size-m::body_width)/2+0.5);}
                for (unsigned int i = 0; i < nodeTagsIn5.size(); i++)   {nodes[nodeTagsIn5[i]].rotate(tau,step,w, m::head_size + m::body_height, m::z_offset+ (m::head_size-m::body_width)/2+1);}
                for (unsigned int i = 0; i < nodeTagsIn3.size(); i++)    {nodes[nodeTagsIn3[i]].rotate(tau,step,-w, m::head_size, m::z_offset+ (m::head_size-m::body_width)/2+0.5);}
                for (unsigned int i = 0; i < nodeTagsIn4.size(); i++)   {nodes[nodeTagsIn4[i]].rotate(tau,step,-w, m::head_size + m::body_height, m::z_offset+ (m::head_size-m::body_width)/2+1);}
        }
        if(step%20 >5 and step%20 <10){
            for (unsigned int i = 0; i < nodeTagsIn2.size(); i++)   {nodes[nodeTagsIn2[i]].rotate(tau,step,-w, m::head_size, m::z_offset+ (m::head_size-m::body_width)/2+0.5);}
                for (unsigned int i = 0; i < nodeTagsIn5.size(); i++)   {nodes[nodeTagsIn5[i]].rotate(tau,step,-w, m::head_size + m::body_height, m::z_offset+ (m::head_size-m::body_width)/2+1);}
                for (unsigned int i = 0; i < nodeTagsIn3.size(); i++)    {nodes[nodeTagsIn3[i]].rotate(tau,step,w, m::head_size, m::z_offset+ (m::head_size-m::body_width)/2+0.5);}
                for (unsigned int i = 0; i < nodeTagsIn4.size(); i++)   {nodes[nodeTagsIn4[i]].rotate(tau,step,w, m::head_size + m::body_height, m::z_offset+ (m::head_size-m::body_width)/2+1);}
        }
        if(step%20 >10 and step%20 <15){
            for (unsigned int i = 0; i < nodeTagsIn2.size(); i++)   {nodes[nodeTagsIn2[i]].rotate(tau,step,-w, m::head_size, m::z_offset+ (m::head_size-m::body_width)/2+0.5);}
                for (unsigned int i = 0; i < nodeTagsIn5.size(); i++)   {nodes[nodeTagsIn5[i]].rotate(tau,step,-w, m::head_size + m::body_height, m::z_offset+ (m::head_size-m::body_width)/2+1);}
                for (unsigned int i = 0; i < nodeTagsIn3.size(); i++)    {nodes[nodeTagsIn3[i]].rotate(tau,step,w, m::head_size, m::z_offset+ (m::head_size-m::body_width)/2+0.5);}
                for (unsigned int i = 0; i < nodeTagsIn4.size(); i++)   {nodes[nodeTagsIn4[i]].rotate(tau,step,w, m::head_size + m::body_height, m::z_offset+ (m::head_size-m::body_width)/2+1);}
        }
        if(step%20 <20 and step%20> 15){
            for (unsigned int i = 0; i < nodeTagsIn2.size(); i++)   {nodes[nodeTagsIn2[i]].rotate(tau,step,w, m::head_size, m::z_offset+ (m::head_size-m::body_width)/2+0.5);}
            for (unsigned int i = 0; i < nodeTagsIn5.size(); i++)   {nodes[nodeTagsIn5[i]].rotate(tau,step,w, m::head_size + m::body_height, m::z_offset+ (m::head_size-m::body_width)/2+1);}
            for (unsigned int i = 0; i < nodeTagsIn3.size(); i++)    {nodes[nodeTagsIn3[i]].rotate(tau,step,-w, m::head_size, m::z_offset+ (m::head_size-m::body_width)/2+0.5);}
            for (unsigned int i = 0; i < nodeTagsIn4.size(); i++)   {nodes[nodeTagsIn4[i]].rotate(tau,step,-w, m::head_size + m::body_height, m::z_offset+ (m::head_size-m::body_width)/2+1);}
    }
        
    }

    // Метод отвечает за запись текущего состояния сетки в снапшот в формате VTK
    void snapshot(unsigned int snap_number) {
        // Сетка в терминах VTK
        vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
        // Точки сетки в терминах VTK
        vtkSmartPointer<vtkPoints> dumpPoints = vtkSmartPointer<vtkPoints>::New();

        // Скалярное поле на точках сетки
        auto smth = vtkSmartPointer<vtkDoubleArray>::New();
        smth->SetName("smth");

        // Векторное поле на точках сетки
        auto vel = vtkSmartPointer<vtkDoubleArray>::New();
        vel->SetName("velocity");
        vel->SetNumberOfComponents(3);

        // Обходим все точки нашей расчётной сетки
        for(unsigned int i = 0; i < nodes.size(); i++) {
            // Вставляем новую точку в сетку VTK-снапшота
            dumpPoints->InsertNextPoint(nodes[i].x, nodes[i].y, nodes[i].z);

            // Добавляем значение векторного поля в этой точке
            double _vel[3] = {nodes[i].vx, nodes[i].vy, nodes[i].vz};
            vel->InsertNextTuple(_vel);

            // И значение скалярного поля тоже
            smth->InsertNextValue(nodes[i].smth);
        }

        // Грузим точки в сетку
        unstructuredGrid->SetPoints(dumpPoints);

        // Присоединяем векторное и скалярное поля к точкам
        unstructuredGrid->GetPointData()->AddArray(vel);
        unstructuredGrid->GetPointData()->AddArray(smth);

        // А теперь пишем, как наши точки объединены в тетраэдры
        for(unsigned int i = 0; i < elements.size(); i++) {
            auto tetra = vtkSmartPointer<vtkTetra>::New();
            tetra->GetPointIds()->SetId( 0, elements[i].nodesIds[0] );
            tetra->GetPointIds()->SetId( 1, elements[i].nodesIds[1] );
            tetra->GetPointIds()->SetId( 2, elements[i].nodesIds[2] );
            tetra->GetPointIds()->SetId( 3, elements[i].nodesIds[3] );
            unstructuredGrid->InsertNextCell(tetra->GetCellType(), tetra->GetPointIds());
        }

        // Создаём снапшот в файле с заданным именем
        string fileName = "tetr3d-step-" + std::to_string(snap_number) + ".vtu";
        vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
        writer->SetFileName(fileName.c_str());
        writer->SetInputData(unstructuredGrid);
        writer->Write();
    }
};


int main(int argc, char **argv)
{
    // Шаг точек по пространству
    double h = 4.0;
    // Шаг по времени

    const unsigned int GMSH_TETR_CODE = 4;

    // Теперь придётся немного упороться:
    // (а) построением сетки средствами gmsh,
    // (б) извлечением данных этой сетки в свой код.

    gmsh::initialize(argc, argv);
  
    // ----------------------- Создание геометрии -----------------------
    gmsh::model::mesh::createGeometry();
    // Голова (Куб)
    gmsh::model::occ::addBox(0, 0, m::z_offset, m::head_size, m::head_size, m::head_size, 1);
  
    // Тело (Куб)
    gmsh::model::occ::addBox((m::head_size - m::body_width) / 2, m::head_size, m::z_offset + (m::head_size-m::body_width)/2, m::body_width, m::body_height, m::body_width, 2);
  
    // Правая рука (Куб)
    gmsh::model::occ::addBox(m::head_size + (m::body_width - m::arm_width) / 2-1.5, m::head_size, m::z_offset+ (m::head_size-m::body_width)/2+0.5, m::arm_width, m::arm_height, m::arm_width, 3);
  
    // Левая рука (Куб)
    gmsh::model::occ::addBox(-(m::body_width + m::arm_width) / 2 + (m::head_size - m::body_width) / 2+0.5, m::head_size, m::z_offset+ (m::head_size-m::body_width)/2+0.5 ,m::arm_width, m::arm_height, m::arm_width, 4);
  
    // Правая нога (Куб)
    gmsh::model::occ::addBox((m::head_size - m::body_width) / 2+ m::leg_width-2, m::head_size + m::body_height, m::z_offset+ (m::head_size-m::body_width)/2+1, m::leg_width, m::leg_height, m::leg_width, 5);
  
    // Левая нога (Куб)
    gmsh::model::occ::addBox((m::head_size + m::body_width) / 2-2, m::head_size + m::body_height, m::z_offset+ (m::head_size-m::body_width)/2+1, m::leg_width, m::leg_height, m::leg_width, 6);

    // ----------------------- Объединение геометрии -----------------------
  
    std::vector<std::pair<int, int>> ov;
    std::vector<std::pair<int, int>> out;
    std::vector<std::vector<std::pair<int, int>>> outDimTags; 
    gmsh::model::occ::fragment({{3,1},{3,2},{3,3},{3,4},{3,5},{3,6}}, ov, out, outDimTags, -1, true, false); // объединение всех тел
  
    gmsh::model::occ::synchronize();


    // Установите размер элементов сетки
    gmsh::option::setNumber("Mesh.CharacteristicLengthMin", m::scale / 2);
    gmsh::option::setNumber("Mesh.CharacteristicLengthMax", m::scale / 1);
  
    // Сгенерируйте сетку (3D)
    gmsh::model::mesh::generate(3);
    //gmsh::write("torus_occ.msh");


    std::vector<double> nodesCoord;
    std::vector<std::size_t> nodeTags;
    std::vector<double> parametricCoord;
    gmsh::model::mesh::getNodes(nodeTags, nodesCoord, parametricCoord);

    // И данные об элементах сетки тоже извлечём, нам среди них нужны только тетраэдры, которыми залит объём



    std::vector<std::size_t>* tetrsNodesTags = nullptr;
    std::vector<int> elementTypes;
    std::vector<std::vector<std::size_t>> elementTags;
    std::vector<std::vector<std::size_t>> elementNodeTags;

    std::vector<std::size_t> nodeTagsIn1;
    std::vector<std::size_t> nodeTagsIn2;
    std::vector<std::size_t> nodeTagsIn3;
    std::vector<std::size_t> nodeTagsIn4;
    std::vector<std::size_t> nodeTagsIn5;
    gmsh::model::mesh::getElements(elementTypes, elementTags, elementNodeTags);

    for(unsigned int i = 0; i < elementTypes.size(); i++) {
        if(elementTypes[i] != GMSH_TETR_CODE)
            continue;
        tetrsNodesTags = &elementNodeTags[i];
    }

    if(tetrsNodesTags == nullptr) {
        cout << "Can not find tetra data. Exiting." << endl;
        gmsh::finalize();
        return -2;
    }

    cout << "The model has " <<  nodeTags.size() << " nodes and " << tetrsNodesTags->size() / 4 << " tetrs." << endl;

    for(unsigned i =0; i<nodeTags.size(); i++){
        if (is_head(nodesCoord, i+1)){nodeTagsIn1.push_back(i);}
        else if (is_body(nodesCoord, i+1)){nodeTagsIn1.push_back(i);}
        else if (is_right_hand(nodesCoord, i+1)){nodeTagsIn2.push_back(i);}
        else if (is_left_hand(nodesCoord, i+1)){nodeTagsIn3.push_back(i);}
        else if (is_right_leg(nodesCoord, i+1)){nodeTagsIn4.push_back(i);}
        else if (is_left_leg(nodesCoord, i+1)){nodeTagsIn5.push_back(i);}
    }

    // На всякий случай проверим, что номера узлов идут подряд и без пробелов
    for(int i = 0; i < nodeTags.size(); ++i) {
        // Индексация в gmsh начинается с 1, а не с нуля. Ну штош, значит так.
        assert(i == nodeTags[i] - 1);
    }
    // И ещё проверим, что в тетраэдрах что-то похожее на правду лежит.
    assert(tetrsNodesTags->size() % 4 == 0);

    // TODO: неплохо бы полноценно данные сетки проверять, да
    CalcMesh mesh(nodesCoord, *tetrsNodesTags);
    double tau = 0.01;
    
    mesh.snapshot(0);
    for(unsigned int step = 1; step < 150; step++) {
        mesh.doTimeStep(tau, step,nodeTagsIn1,nodeTagsIn2,nodeTagsIn3,nodeTagsIn4,nodeTagsIn5);
        mesh.snapshot(step);
    }
    gmsh::finalize();
    return 0;
}




/*    int pg1 = gmsh::model::addPhysicalGroup(3, {1,2}); // head + body
    gmsh::model::setPhysicalName(3, pg1, "HB"); 

    int pg2 = gmsh::model::addPhysicalGroup(3, {3}); 
    gmsh::model::setPhysicalName(3, pg2, "RH"); 

    int pg3 = gmsh::model::addPhysicalGroup(3, {4});
    gmsh::model::setPhysicalName(3, pg3, "LH"); 

    int pg4 = gmsh::model::addPhysicalGroup(3, {5}); 
    gmsh::model::setPhysicalName(3, pg4, "RL"); 

    int pg5 =  gmsh::model::addPhysicalGroup(3, {6});
    gmsh::model::setPhysicalName(3, pg5, "LL"); 
    std::cout << pg3;
    // ----------------------- Меширование -----------------------
  */

  

    /*std::vector<int> elementTypes1;
    std::vector<std::vector<std::size_t>> elementTags1;
    std::vector<std::vector<std::size_t>> elementNodeTags1;

    std::vector<int> elementTypes2;
    std::vector<std::vector<std::size_t>> elementTags2;
    std::vector<std::vector<std::size_t>> elementNodeTags2;

    std::vector<int> elementTypes3;
    std::vector<std::vector<std::size_t>> elementTags3;
    std::vector<std::vector<std::size_t>> elementNodeTags3;

    std::vector<int> elementTypes4;
    std::vector<std::vector<std::size_t>> elementTags4;
    std::vector<std::vector<std::size_t>> elementNodeTags4;

    std::vector<int> elementTypes5;
    std::vector<std::vector<std::size_t>> elementTags5;
    std::vector<std::vector<std::size_t>> elementNodeTags5;

    gmsh::model::mesh::getElements(elementTypes1, elementTags1, elementNodeTags1, 1);
    gmsh::model::mesh::getElements(elementTypes2, elementTags2, elementNodeTags2, 2);
    gmsh::model::mesh::getElements(elementTypes3, elementTags3, elementNodeTags3, 3);
    gmsh::model::mesh::getElements(elementTypes4, elementTags4, elementNodeTags4, 4);
    gmsh::model::mesh::getElements(elementTypes5, elementTags5, elementNodeTags5, 5);*/

    /*std::vector<std::size_t> nodeTagsIn1;
    std::vector<double> nodeCoordinatesIn1;
    gmsh::model::mesh::getNodes(nodeTagsIn1, nodeCoordinatesIn1, pg1);

    std::vector<std::size_t> nodeTagsIn5;
    std::vector<double> nodeCoordinatesIn5;
    gmsh::model::mesh::getNodes(nodeTagsIn5, nodeCoordinatesIn5, pg5);

    std::vector<std::size_t> nodeTagsIn2;
    std::vector<double> nodeCoordinatesIn2;
    gmsh::model::mesh::getNodes(nodeTagsIn2, nodeCoordinatesIn2, pg2);

    std::vector<std::size_t> nodeTagsIn3;
    std::vector<double> nodeCoordinatesIn3;
    gmsh::model::mesh::getNodes(nodeTagsIn3, nodeCoordinatesIn3, pg3);

    std::vector<std::size_t> nodeTagsIn4;
    std::vector<double> nodeCoordinatesIn4;
    gmsh::model::mesh::getNodes(nodeTagsIn4, nodeCoordinatesIn4, pg4);*/

    //std::array<int, 5> physicalGroups = {pg1, pg2, pg3, pg4, pg5};

    //