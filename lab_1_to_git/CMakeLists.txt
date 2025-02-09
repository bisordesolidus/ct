cmake_minimum_required(VERSION 3.11)

set(CMAKE_CXX_FLAGS "-std=c++11 -g")

project(gmshdemo)
add_executable(tor tor.cpp)
add_executable(poor_attempts_to_task2 poor_attempts_to_task2.cpp)
add_executable(task2_correct task2_correct.cpp)
# --- Подключаем gmsh, тут придется перебирать варианты 

# - если вам повезло, то сработает простой путь
#include_directories(${GMSH_INCLUDE_DIR})

# - если у вас Ubuntu, то возможно, что приехал немного кривой gmsh из репозитория
# - при этом достаточно указать системные пути, куда система положила нужные файлы
# - закомментируйте строку 14 и раскомментируйте эти шесть:
add_library(gmsh SHARED IMPORTED)
set_target_properties(gmsh PROPERTIES
    IMPORTED_LOCATION "/usr/lib/x86_64-linux-gnu/libgmsh.so"
    INTERFACE_INCLUDE_DIRECTORIES "/usr/include"
)
set(GMSH_LIBRARIES gmsh)

# - если предыдущий вариант не сработал, собирайте gmsh из исходников
# - после этого в предыдущем варианте заменяйте системные пути на те, что вы сделали сами

# --- Закончили с gmsh, остальное править не обязательно (не поможет)

# Укажем библиотеки для линковки
target_link_libraries(tor ${GMSH_LIBRARIES})
target_link_libraries(poor_attempts_to_task2 ${GMSH_LIBRARIES})
target_link_libraries(task2_correct ${GMSH_LIBRARIES})