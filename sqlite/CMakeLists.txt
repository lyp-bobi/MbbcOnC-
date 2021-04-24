
add_executable(shell_sqlite sqlite3.c shell.c)

set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DSQLITE_THREADSAFE=2 -DSQLITE_OMIT_LOAD_EXTENSION -lm -pthread")
set (CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -DSQLITE_THREADSAFE=2 -DSQLITE_OMIT_LOAD_EXTENSION -lm -pthread")

#if(WIN32)
#    add_library(sqlite3_lib SHARED sqlite3.dll)
#    set_target_properties(sqlite3_lib PROPERTIES LINKER_LANGUAGE CXX)
#else()
add_library(sqlite3_lib SHARED sqlite3.c shell.c)
#endif(WIN32)

install(TARGETS sqlite3_lib
        RUNTIME DESTINATION ${SIDX_BIN_DIR}
        LIBRARY DESTINATION ${SIDX_LIB_DIR}
        ARCHIVE DESTINATION ${SIDX_LIB_DIR})
