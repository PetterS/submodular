# Author: petter.strandmark@gmail.com (Petter Strandmark)


ADD_LIBRARY(QPBO QPBO.cpp QPBO_extra.cpp QPBO_maxflow.cpp QPBO_postprocessing.cpp)

INSTALL(TARGETS QPBO
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)
