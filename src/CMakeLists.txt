CMAKE_MINIMUM_REQUIRED ( VERSION 2.8.0 )
PROJECT( NBTISym )

FILE(TO_CMAKE_PATH "${CMAKE_CURRENT_SOURCE_DIR}" ADD_CURRENT_CMAKE_MODULE_PATH)

list(APPEND CMAKE_MODULE_PATH "${FLEX_PATH}")

find_package(FLEX REQUIRED)
find_package(BISON REQUIRED)

SET( SRC_TRACKING
  easylogging++.h
  ezOptionParser.hpp
  FlexLexer.h
  main.cpp
  ProcessorState.hpp
  utils.hpp
  parser/parser.yy
  parser/scanner.ll
  parser/driver.cc
  parser/driver.hh
)

SOURCE_GROUP( "Parser files" FILES 
  parser/parser.yy
  parser/scanner.ll 
  parser/driver.cc
  parser/driver.hh
)

if(WIN32)
	SET( WIN_COMPAT_FLAG --wincompat )
endif(WIN32)

BISON_TARGET(Parser parser/parser.yy ${CMAKE_CURRENT_SOURCE_DIR}/parser/parser.cc)
FLEX_TARGET(Scanner parser/scanner.ll  ${CMAKE_CURRENT_SOURCE_DIR}/parser/scanner.cc COMPILE_FLAGS ${WIN_COMPAT_FLAG})
ADD_FLEX_BISON_DEPENDENCY(Scanner Parser)

add_executable( nbtisym ${SRC_TRACKING} ${BISON_Parser_OUTPUTS} ${FLEX_Scanner_OUTPUTS} )

TARGET_LINK_LIBRARIES( ${CMAKE_PROJECT_NAME} )