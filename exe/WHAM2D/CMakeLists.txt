include_directories(.)

add_executable (WHAM2D_makeDOS makeDOS.cpp)
target_link_libraries(WHAM2D_makeDOS ${Boost_LIBRARIES})

add_executable (WHAM2D_makeDOSJK makeDOSJK.cpp)
target_link_libraries(WHAM2D_makeDOSJK ${Boost_LIBRARIES})

add_executable (WHAM2D_reweight_constK reweight_constK.cpp)
target_link_libraries(WHAM2D_reweight_constK ${Boost_LIBRARIES})

add_executable (WHAM2D_reweightdT_constK reweightdT_constK.cpp)
target_link_libraries(WHAM2D_reweightdT_constK ${Boost_LIBRARIES})

add_executable (WHAM2D_logWeights logWeights.cpp)
target_link_libraries(WHAM2D_logWeights ${Boost_LIBRARIES})
