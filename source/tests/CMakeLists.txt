# Author: petter.strandmark@gmail.com (Petter Strandmark)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/thirdparty/Catch)

MACRO (SUBMODULAR_TEST NAME)
	ADD_EXECUTABLE(${NAME} 
	               ${NAME}.cpp
	               ${SPII_HEADERS})

	TARGET_LINK_LIBRARIES(${NAME} pseudoboolean)
	TARGET_LINK_LIBRARIES(${NAME} ${SUBMODULAR_LIBRARY_DEPENDENCIES})
	SET_PROPERTY(TARGET ${NAME}         PROPERTY FOLDER "Tests")
	ADD_TEST(NAME ${NAME}
	         COMMAND ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${NAME}
	         WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
ENDMACRO ()

FILE(GLOB TEST_FILES test_*.cpp)
FOREACH (TEST_FILE ${TEST_FILES})
	GET_FILENAME_COMPONENT(TEST_NAME ${TEST_FILE} NAME_WE)
	MESSAGE("-- Adding test: " ${TEST_NAME})
	SUBMODULAR_TEST(${TEST_NAME})
ENDFOREACH()

#
# Valgrind (memcheck) testing.
#
find_program(VALGRIND NAMES valgrind PATH /usr/bin /usr/local/bin)
if (VALGRIND)
	message("-- Valgrind found; configuring memcheck tests.")
	macro(VALGRIND_TEST NAME)
		add_test(NAME memcheck_${NAME}
		         COMMAND ${VALGRIND} --error-exitcode=1 --leak-check=full ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_${NAME})
	endmacro(VALGRIND_TEST)

	# Run shorter tests though memcheck.
	VALGRIND_TEST(function)
	VALGRIND_TEST(function_serializer)
	VALGRIND_TEST(solver)
	VALGRIND_TEST(spii)
	VALGRIND_TEST(term)
	VALGRIND_TEST(term_factory)
endif (VALGRIND)

#
# Cppcheck testing.
#
find_program(CPPCHECK NAMES cppcheck PATH /usr/bin /usr/local/bin)
if (CPPCHECK)
	message("-- Cppcheck found; configuring test.")
	add_test(NAME cppcheck
	         COMMAND ${CPPCHECK} ${CMAKE_SOURCE_DIR} -i ${CMAKE_SOURCE_DIR}/thirdparty/ -q --error-exitcode=1)
endif (CPPCHECK)
