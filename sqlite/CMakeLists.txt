add_executable(sqlite sqlite3.c shell.c)

add_library(${SIDX_LIB_NAME} SHARED sqlite3.c shell.c)