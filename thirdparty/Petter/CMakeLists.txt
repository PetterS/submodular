# Author: petter.strandmark@gmail.com (Petter Strandmark)


ADD_LIBRARY(petter Petter-Color.cc)

INSTALL(TARGETS petter
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)
