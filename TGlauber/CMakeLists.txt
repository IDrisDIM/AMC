project(TGlauberLib)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(Source
	TGlauNucleon.cc
	TGlauNucleus.cc
	TGlauberMC.cc)

add_library(TGlauberLib
	${Source})
target_link_libraries(TGlauberLib
	${ROOT_LIBRARIES})
