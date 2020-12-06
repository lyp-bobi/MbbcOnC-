include_directories(../include)

find_package(Threads)
set (DIR geometry)
set (SOURCES 
        Intersection
        )

foreach (test ${SOURCES})
    add_executable(test-${DIR}-${test} ${DIR}/${test}.cc)
    target_link_libraries(test-${DIR}-${test} ${SIDX_LIB_NAME})
    target_link_libraries(test-${DIR}-${test} ${CMAKE_THREAD_LIBS_INIT})
endforeach()
#
#set (DIR rtree)
#set (SOURCES
#        Exhaustive
#        Generator
#        RTreeBulkLoad
#        RTreeLoad
#        RTreeQuery)
#
#
#foreach (test ${SOURCES})
#    add_executable(test-${DIR}-${test} ${DIR}/${test}.cc)
#    target_link_libraries(test-${DIR}-${test} ${SIDX_LIB_NAME})
#endforeach()
#
#set (DIR mvrtree)
#set (SOURCES
#        Exhaustive
#        Generator
#        MVRTreeLoad
#        MVRTreeQuery
#        )
#
#
#foreach (test ${SOURCES})
#    add_executable(test-${DIR}-${test} ${DIR}/${test}.cc)
#    target_link_libraries(test-${DIR}-${test} ${SIDX_LIB_NAME})
#endforeach()
#
#set (DIR tprtree)
#set (SOURCES
#        Exhaustive
#        Generator
#        TPRTreeLoad
#        TPRTreeQuery
#        )
#
#foreach (test ${SOURCES})
#    add_executable(test-${DIR}-${test} ${DIR}/${test}.cc ${DIR}/RandomGenerator.cc)
#    target_link_libraries(test-${DIR}-${test} ${SIDX_LIB_NAME})
#endforeach()


set (DIR mbc)
set (SOURCES
        #mbcTest
        leafsize
        knnTest
        densityTest
        TrajStoreTest
        optimizedTest
        tb3drTest
        e2eTest
        e2eTestGT
        splitTest
        radiusEstimator
        splitQueryTest
        rangeTest
        simpliTest
        speedEstimator
        changeKTest
        e2eRangeTest
        indexTest
        dump
        secondaryTest
        )

foreach (test ${SOURCES})
    add_executable(test-${DIR}-${test} ${DIR}/${test}.cc)
    target_link_libraries(test-${DIR}-${test} ${SIDX_LIB_NAME})
    target_link_libraries(test-${DIR}-${test} ${CMAKE_THREAD_LIBS_INIT})
endforeach()

set (DIR testx)
set (SOURCES
        regressTest
        knnTest
        #knn-mt
        dump
        )

foreach (test ${SOURCES})
    add_executable(test-${DIR}-${test} ${DIR}/${test}.cc testx/knn-mt.cc)
    target_link_libraries(test-${DIR}-${test} ${SIDX_LIB_NAME})
    target_link_libraries(test-${DIR}-${test} ${CMAKE_THREAD_LIBS_INIT})
endforeach()