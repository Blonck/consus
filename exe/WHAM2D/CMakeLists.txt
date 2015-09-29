include_directories(.)

add_executable (WHAM2D_makeDOSandMicroMean makeDOSandMicroMean.cpp)
target_link_libraries(WHAM2D_makeDOSandMicroMean ${Boost_LIBRARIES})
