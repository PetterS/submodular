# Author: petter.strandmark@gmail.com (Petter Strandmark)

configure_file(generators/generators.txt ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/generators/generators.txt COPYONLY)
configure_file(generators/newgenerators.txt ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/generators/newgenerators.txt COPYONLY)

configure_file(batchrun3.cmd ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/batchrun3.cmd COPYONLY)
configure_file(batchrun4.cmd ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/batchrun4.cmd COPYONLY)
configure_file(check_correctness.cmd ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/check_correctness.cmd COPYONLY)
configure_file(check_correctness.cmd ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/check_correctness.cmd COPYONLY)
configure_file("examples from paper.cmd" "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/examples from paper.cmd" COPYONLY)


configure_file(parse_batchrun.py ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/parse_batchrun.py COPYONLY)
configure_file(plot_batchrun.m ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/plot_batchrun.m COPYONLY)
configure_file(plot_batchrun_new.m ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/plot_batchrun_new.m COPYONLY)
