include_directories(.)

add_executable (reweightE reweightE.cpp)
target_link_libraries(reweightE ${Boost_LIBRARIES})

add_executable (NVT2D_makeDOS makeDOS.cpp)
target_link_libraries(NVT2D_makeDOS ${Boost_LIBRARIES})

add_executable (NVT2D_makeHists makeHists.cpp)
target_link_libraries(NVT2D_makeHists ${Boost_LIBRARIES})

add_executable (NVT2D_makeOverlap makeOverlap.cpp)
target_link_libraries(NVT2D_makeOverlap ${Boost_LIBRARIES})

add_executable (NVT2D_makeMicroMeans makeMicroMeans.cpp)
target_link_libraries(NVT2D_makeMicroMeans ${Boost_LIBRARIES})

add_executable (NVT2D_makeDOSJK makeDOSJK.cpp)
target_link_libraries(NVT2D_makeDOSJK ${Boost_LIBRARIES})

add_executable (NVT2D_reweight_constK reweight_constK.cpp)
target_link_libraries(NVT2D_reweight_constK ${Boost_LIBRARIES})

add_executable (NVT2D_reweightdT_constK reweightdT_constK.cpp)
target_link_libraries(NVT2D_reweightdT_constK ${Boost_LIBRARIES})

add_executable (NVT2D_reweightdT2_constK reweightdT2_constK.cpp)
target_link_libraries(NVT2D_reweightdT2_constK ${Boost_LIBRARIES})

add_executable (NVT2D_reweight_contour reweight_contour.cpp)
target_link_libraries(NVT2D_reweight_contour ${Boost_LIBRARIES})

add_executable (NVT2D_reweightdT_contour reweightdT_contour.cpp)
target_link_libraries(NVT2D_reweightdT_contour ${Boost_LIBRARIES})

add_executable (NVT2D_searchMinMax searchMinMax.cpp)
target_link_libraries(NVT2D_searchMinMax ${Boost_LIBRARIES})
