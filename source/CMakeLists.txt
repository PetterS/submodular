# Author: petter.strandmark@gmail.com (Petter Strandmark)

ADD_EXECUTABLE(submodular
               main_program.cpp
               submodular_tests.cpp
               ${SUBMODULAR_HEADERS})
TARGET_LINK_LIBRARIES(submodular pseudoboolean)
TARGET_LINK_LIBRARIES(submodular ${SUBMODULAR_LIBRARY_DEPENDENCIES})


ADD_TEST(NAME library COMMAND submodular WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
