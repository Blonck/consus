include_directories(.)

add_executable (WHAM2D_makeDOSandMicroMean makeDOSandMicroMean.cpp)
target_link_libraries(WHAM2D_makeDOSandMicroMean ${Boost_LIBRARIES})

add_executable (WHAM2D_reweight_constK reweight_constK.cpp)
target_link_libraries(WHAM2D_reweight_constK ${Boost_LIBRARIES})

add_executable (WHAM2D_logWeights logWeights.cpp)
target_link_libraries(WHAM2D_logWeights ${Boost_LIBRARIES})
