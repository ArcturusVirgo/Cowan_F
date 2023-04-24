# The following CMake script was shamelessly stolen from the CMake mailing list
#      https://cmake.org/pipermail/cmake/2009-July/030788.html
# It is used to compare the output of an test executable against a known-good
# ("blessed") output

# some argument checking:
# test_cmd is the command to run with all its arguments
IF(NOT test_cmd)
   MESSAGE(FATAL_ERROR "Variable test_cmd not defined")
ENDIF(NOT test_cmd)

# output_blessed contains the name of the "blessed" output file
IF(NOT output_blessed)
   MESSAGE(FATAL_ERROR "Variable output_blessed not defined")
ENDIF(NOT output_blessed)

# output_test contains the name of the output file the test_cmd will produce
IF(NOT output_test)
   MESSAGE(FATAL_ERROR "Variable output_test not defined")
ENDIF(NOT output_test)

# convert the space-separated string to a list
SEPARATE_ARGUMENTS(test_args)
MESSAGE(${test_args})

EXECUTE_PROCESS(
   COMMAND ${test_cmd} ${test_args}
)

EXECUTE_PROCESS(
   COMMAND ${CMAKE_COMMAND} -E compare_files ${output_blessed} ${output_test}
   RESULT_VARIABLE test_not_successful
)

IF(test_not_successful)
   MESSAGE("Error message: ${test_not_successful}")
   MESSAGE(SEND_ERROR "${output_test} does not match ${output_blessed}!")
ENDIF(test_not_successful)
