cmake_minimum_required(VERSION 3.16)
project(hq_hadronizer_cpp LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(${PROJECT_NAME}
  src/main.cpp
  src/RNG.cpp
  src/Lorentz.cpp
  src/PhysicsLists.cpp
  src/Hadronizer.cpp
  src/IO.cpp
  src/Recomb2.cpp
  src/Recomb3.cpp
  src/RecombinationTable.cpp
  src/Frag.cpp
)

target_include_directories(${PROJECT_NAME} PRIVATE include)

# Optional: warnings
if (MSVC)
  target_compile_options(${PROJECT_NAME} PRIVATE /W4)
else()
  target_compile_options(${PROJECT_NAME} PRIVATE -Wall -Wextra -Wpedantic)
endif()

