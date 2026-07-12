cmake_minimum_required(VERSION 3.16)

project(LinearAlgebra LANGUAGES CXX)

add_executable(LinearAlgebra
    LinearAlgebra/Matrix.cpp
    LinearAlgebra/Source.cpp
)

target_compile_features(LinearAlgebra PRIVATE cxx_std_17)
target_include_directories(LinearAlgebra PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/LinearAlgebra
)

if (MSVC)
    target_compile_options(LinearAlgebra PRIVATE /W4)
else()
    target_compile_options(LinearAlgebra PRIVATE -Wall -Wextra -pedantic)
endif()
