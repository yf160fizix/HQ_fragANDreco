cmake_minimum_required(VERSION 3.16)
project(hq_hadronizer_cpp LANGUAGES CXX)

find_package(Boost REQUIRED COMPONENTS program_options)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(ISRM_STRICT_WARNINGS "Enable additional conversion and shadow warnings" OFF)
option(ISRM_ENABLE_SANITIZERS "Enable AddressSanitizer and UndefinedBehaviorSanitizer" OFF)

include(CheckCXXSourceCompiles)
check_cxx_source_compiles(
  "#include <cmath>\nint main() { return std::cyl_bessel_k(2.0, 1.0) > 0.0 ? 0 : 1; }"
  ISRM_HAS_STD_CYL_BESSEL_K
)
if(ISRM_HAS_STD_CYL_BESSEL_K)
  message(STATUS "Using std::cyl_bessel_k")
else()
  message(STATUS "std::cyl_bessel_k unavailable; using Boost.Math fallback")
endif()

add_library(hadronization STATIC
  src/RNG.cpp
  src/Lorentz.cpp
  src/Hadronizer.cpp
  src/IO.cpp
  src/Recomb2.cpp
  src/Recomb3.cpp
  src/RecombKinematics.cpp
  src/RecombChannel.cpp
  src/RecombTable.cpp
  src/Frag.cpp
  src/CharmFeeddown.cpp
  src/CharmStateTable.cpp
)
target_include_directories(hadronization PUBLIC include)
target_include_directories(hadronization SYSTEM PRIVATE ${Boost_INCLUDE_DIRS})
if(ISRM_HAS_STD_CYL_BESSEL_K)
  target_compile_definitions(hadronization PRIVATE ISRM_HAS_STD_CYL_BESSEL_K=1)
endif()

add_executable(${PROJECT_NAME} src/main.cpp)
target_link_libraries(${PROJECT_NAME} PRIVATE hadronization Boost::program_options)

if(MSVC)
  target_compile_options(hadronization PRIVATE /W4)
  target_compile_options(${PROJECT_NAME} PRIVATE /W4)
else()
  set(ISRM_WARNING_FLAGS -Wall -Wextra -Wpedantic)
  if(ISRM_STRICT_WARNINGS)
    list(APPEND ISRM_WARNING_FLAGS -Wshadow -Wconversion -Wsign-conversion)
  endif()
  target_compile_options(hadronization PRIVATE ${ISRM_WARNING_FLAGS})
  target_compile_options(${PROJECT_NAME} PRIVATE ${ISRM_WARNING_FLAGS})

  if(ISRM_ENABLE_SANITIZERS)
    set(ISRM_SANITIZER_FLAGS -fsanitize=address,undefined -fno-omit-frame-pointer)
    target_compile_options(hadronization PRIVATE ${ISRM_SANITIZER_FLAGS})
    target_compile_options(${PROJECT_NAME} PRIVATE ${ISRM_SANITIZER_FLAGS})
    target_link_options(hadronization INTERFACE -fsanitize=address,undefined)
  endif()
endif()

include(CTest)
if(BUILD_TESTING)
  add_executable(test_reco_feeddown tests/test_reco_feeddown.cpp)
  target_link_libraries(test_reco_feeddown PRIVATE hadronization)
  add_test(NAME recombination_feeddown COMMAND test_reco_feeddown)

  add_executable(test_lorentz tests/test_lorentz.cpp)
  target_link_libraries(test_lorentz PRIVATE hadronization)
  add_test(NAME lorentz_round_trip COMMAND test_lorentz)

  add_executable(test_io tests/test_io.cpp)
  target_link_libraries(test_io PRIVATE hadronization)
  add_test(
    NAME oscar_io
    COMMAND test_io
            ${CMAKE_CURRENT_SOURCE_DIR}/tests/data/standard_input.oscar
  )

  add_executable(test_origin tests/test_origin.cpp)
  target_link_libraries(test_origin PRIVATE hadronization)
  add_test(
    NAME hadronization_origin
    COMMAND test_origin
            ${CMAKE_CURRENT_SOURCE_DIR}/data/charm_meson_states.csv
            ${CMAKE_CURRENT_SOURCE_DIR}/data/charm_baryon_states.csv
            ${CMAKE_CURRENT_SOURCE_DIR}/data/recomb_c_raw_M020_B0267.dat
            ${CMAKE_CURRENT_SOURCE_DIR}/data/max_wigner_c_M020_B0267.dat
            ${CMAKE_CURRENT_SOURCE_DIR}/tests/data/regression_input.oscar
  )

  add_executable(test_peterson_parameters tests/test_peterson_parameters.cpp)
  target_link_libraries(test_peterson_parameters PRIVATE hadronization)
  add_test(NAME peterson_parameters COMMAND test_peterson_parameters)

  add_executable(test_gamma_hb tests/test_gamma_hb.cpp)
  target_link_libraries(test_gamma_hb PRIVATE hadronization)
  add_test(
    NAME heavy_baryon_suppression
    COMMAND test_gamma_hb
            ${CMAKE_CURRENT_SOURCE_DIR}/data/charm_meson_states.csv
            ${CMAKE_CURRENT_SOURCE_DIR}/data/charm_baryon_states.csv
            ${CMAKE_CURRENT_SOURCE_DIR}/data/recomb_c_raw_M020_B0267.dat
            ${CMAKE_CURRENT_SOURCE_DIR}/data/max_wigner_c_M020_B0267.dat
  )

  add_executable(test_tchem tests/test_tchem.cpp)
  target_link_libraries(test_tchem PRIVATE hadronization)
  add_test(
    NAME tchem_chemistry
    COMMAND test_tchem
            ${CMAKE_CURRENT_SOURCE_DIR}/data/charm_meson_states.csv
            ${CMAKE_CURRENT_SOURCE_DIR}/data/charm_baryon_states.csv
  )

  if(ISRM_HAS_STD_CYL_BESSEL_K)
    add_executable(test_bessel_chemistry tests/test_bessel_chemistry.cpp)
    target_link_libraries(test_bessel_chemistry PRIVATE hadronization)
    target_include_directories(test_bessel_chemistry SYSTEM PRIVATE
                               ${Boost_INCLUDE_DIRS})
    add_test(
      NAME bessel_chemistry_agreement
      COMMAND test_bessel_chemistry
              ${CMAKE_CURRENT_SOURCE_DIR}/data/charm_meson_states.csv
              ${CMAKE_CURRENT_SOURCE_DIR}/data/charm_baryon_states.csv
    )
  endif()

  if(MSVC)
    target_compile_options(test_reco_feeddown PRIVATE /W4)
    target_compile_options(test_lorentz PRIVATE /W4)
    target_compile_options(test_io PRIVATE /W4)
    target_compile_options(test_origin PRIVATE /W4)
    target_compile_options(test_peterson_parameters PRIVATE /W4)
    target_compile_options(test_gamma_hb PRIVATE /W4)
    target_compile_options(test_tchem PRIVATE /W4)
    if(TARGET test_bessel_chemistry)
      target_compile_options(test_bessel_chemistry PRIVATE /W4)
    endif()
  else()
    target_compile_options(test_reco_feeddown PRIVATE ${ISRM_WARNING_FLAGS})
    target_compile_options(test_lorentz PRIVATE ${ISRM_WARNING_FLAGS})
    target_compile_options(test_io PRIVATE ${ISRM_WARNING_FLAGS})
    target_compile_options(test_origin PRIVATE ${ISRM_WARNING_FLAGS})
    target_compile_options(test_peterson_parameters PRIVATE ${ISRM_WARNING_FLAGS})
    target_compile_options(test_gamma_hb PRIVATE ${ISRM_WARNING_FLAGS})
    target_compile_options(test_tchem PRIVATE ${ISRM_WARNING_FLAGS})
    if(TARGET test_bessel_chemistry)
      target_compile_options(test_bessel_chemistry PRIVATE ${ISRM_WARNING_FLAGS})
    endif()
  endif()

  configure_file(tests/run_regression.cmake.in
                 ${CMAKE_CURRENT_BINARY_DIR}/run_regression.cmake @ONLY)
  add_test(NAME fixed_seed_regression
           COMMAND ${CMAKE_COMMAND} -P
                   ${CMAKE_CURRENT_BINARY_DIR}/run_regression.cmake)

  configure_file(tests/run_cli_tests.cmake.in
                 ${CMAKE_CURRENT_BINARY_DIR}/run_cli_tests.cmake @ONLY)
  add_test(NAME cli_mode_and_tchem
           COMMAND ${CMAKE_COMMAND} -P
                   ${CMAKE_CURRENT_BINARY_DIR}/run_cli_tests.cmake)
endif()
