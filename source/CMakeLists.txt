# Author: petter.strandmark@gmail.com (Petter Strandmark)

ADD_EXECUTABLE(submodular
               main_program.cpp
               ${SUBMODULAR_HEADERS})
TARGET_LINK_LIBRARIES(submodular pseudoboolean)
TARGET_LINK_LIBRARIES(submodular ${SUBMODULAR_LIBRARY_DEPENDENCIES})
