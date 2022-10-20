add_subdirectory(SPSDict)

add_executable(WiggleKiller)

target_include_directories(WiggleKiller PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/.
    ${CMAKE_CURRENT_SOURCE_DIR}/SPSDict
    SYSTEM ${ROOT_INCLUDE_DIRS}
)

target_sources(WiggleKiller PRIVATE
    CubicSpline.h
    CubicSpline.cpp
    CutHandler.h
    CutHandler.cpp
    FP_kinematics.h
    FP_kinematics.cpp
    MassLookup.h
    MassLookup.cpp
    WiggleCorrector.h
    WiggleCorrector.cpp
    WiggleKiller.h
    WiggleKiller.cpp
    main.cpp
)

target_link_libraries(WiggleKiller SPSDict ${ROOT_LIBRARIES})
set_target_properties(WiggleKiller PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${WK_BINARY_DIR})