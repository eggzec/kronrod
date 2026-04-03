cmake_minimum_required(VERSION 3.10)

project(kronrod
    VERSION 0.1.0
    LANGUAGES C
)

# ----------------------------
# Library sources
# ----------------------------

set(KRONROD_SRC
    src/kronrod.c
)

# ----------------------------
# Static library
# ----------------------------

add_library(kronrod STATIC ${KRONROD_SRC})

target_include_directories(kronrod
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
        $<INSTALL_INTERFACE:include>
)

set_target_properties(kronrod PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION 0
)

# ----------------------------
# Shared library
# ----------------------------

add_library(kronrod_shared SHARED ${KRONROD_SRC})

target_include_directories(kronrod_shared
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
        $<INSTALL_INTERFACE:include>
)

set_target_properties(kronrod_shared PROPERTIES
    OUTPUT_NAME kronrod
    VERSION ${PROJECT_VERSION}
    SOVERSION 0
)

# ----------------------------
# Installation
# ----------------------------

install(
    FILES src/kronrod.h
    DESTINATION include
)

install(
    TARGETS kronrod kronrod_shared
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
)

# ----------------------------
# Tests
# ----------------------------

option(KRONROD_BUILD_TESTS "Build C tests" ON)

if(KRONROD_BUILD_TESTS)

    enable_testing()

    add_executable(test_kronrod
        tests/test_kronrod.c
    )

    target_link_libraries(test_kronrod
        kronrod
        $<$<NOT:$<PLATFORM_ID:Windows>>:m>
    )

    add_test(
        NAME kronrod_tests
        COMMAND test_kronrod
    )

endif()

message(STATUS "KRONROD C library CMake setup complete.")
