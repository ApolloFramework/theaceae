
add_library(thcCFDlib thcCfdInput.cpp thcFunctions.cpp thcMesh.cpp)
target_link_libraries(thcCFDlib Camellia)

add_executable(thcCFD thcCFD.cpp)
target_link_libraries(thcCFD thcCFDlib)

set (THCCFD_HEADERS thcCfdInput.h thcMesh.h thcFunctions.h)
