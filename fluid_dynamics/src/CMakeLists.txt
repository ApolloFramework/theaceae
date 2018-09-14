
add_library(thcCFDlib thcCfdInput.cpp thcUtils.cpp thcMesh.cpp)
target_link_libraries(thcCFDlib Camellia)

add_executable(thcCFD thcCFD.cpp)
target_link_libraries(thcCFD thcCFDlib)
