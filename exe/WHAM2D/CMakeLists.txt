add_subdirectory(NVT)

add_executable (WHAM2D_logWeights logWeights.cpp)
target_link_libraries(WHAM2D_logWeights ${Boost_LIBRARIES})
