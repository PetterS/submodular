# Author: petter.strandmark@gmail.com (Petter Strandmark)


ADD_LIBRARY(maxflow3 graph.cpp maxflow.cpp)

INSTALL(TARGETS maxflow3
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)
