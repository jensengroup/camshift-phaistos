# Add directories to the list of includes
include_directories (${phaistos_SOURCE_DIR}/src ${phaistos_SOURCE_DIR}/src/models/backbone_dbn ${CMAKE_CURRENT_SOURCE_DIR}/../src ${Boost_INCLUDE_DIRS}) 

# Output binaries to global bin dir
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/test)

# Create executable
add_executable(test_camshift test_camshift.cpp)

# Link executable with library
target_link_libraries(test_camshift libphaistos dbn ${LAPACK_LIBRARY} ${BLAS_LIBRARY})
