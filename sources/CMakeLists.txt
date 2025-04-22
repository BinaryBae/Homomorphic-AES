set(TFHE_PREFIX /home/rudy/libraries/tfhe)
set(APP_NAME Homo-AES)
set(SOME_NAME Homo-AES)
set(MAIN1 ${CMAKE_CURRENT_LIST_DIR}/main.cpp)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

set(APP_SRCS
  ${CMAKE_CURRENT_LIST_DIR}/homo_aes.cpp
  ${CMAKE_CURRENT_LIST_DIR}/keyswitch.cpp
  ${CMAKE_CURRENT_LIST_DIR}/tables.cpp
  ${CMAKE_CURRENT_LIST_DIR}/utils.cpp
  ${CMAKE_CURRENT_LIST_DIR}/tree_mvb.cpp
  )

set(APP_HEADERS
  ${CMAKE_CURRENT_LIST_DIR}/homo_aes.h
  ${CMAKE_CURRENT_LIST_DIR}/keyswitch.h
  ${CMAKE_CURRENT_LIST_DIR}/tables.h
  ${CMAKE_CURRENT_LIST_DIR}/utils.h
  ${CMAKE_CURRENT_LIST_DIR}/tree_mvb.h

  )

add_executable(${APP_NAME} ${APP_SRCS} ${MAIN1} ${APP_HEADERS})

include_directories(${TFHE_PREFIX}/src/include)

find_package(OpenMP)

target_link_libraries(${APP_NAME} ${TFHE_PREFIX}/build/libtfhe/libtfhe-spqlios-avx.so OpenMP::OpenMP_CXX)
target_link_libraries(${SOME_NAME} ${TFHE_PREFIX}/build/libtfhe/libtfhe-spqlios-avx.so OpenMP::OpenMP_CXX)

set_target_properties(${APP_NAME}
  PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
  )

