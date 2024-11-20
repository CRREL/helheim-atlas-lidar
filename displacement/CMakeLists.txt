cmake_minimum_required(VERSION 3.16)

project(ATLAS
    VERSION 1.0
    LANGUAGES CXX
)

find_package(PDAL 2.7 REQUIRED CONFIG)
find_package(GDAL REQUIRED CONFIG)
find_package(Eigen3 REQUIRED CONFIG)

add_executable(atlas
    App.cpp
    Atlas.cpp
    Draw.cpp
    Grid.cpp
)
target_link_libraries(atlas
    PRIVATE
        pdalcpp
        Eigen3::Eigen
        GDAL::GDAL
)

target_compile_features(atlas
    PRIVATE
        cxx_std_11
)
