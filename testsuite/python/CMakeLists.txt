#
# Copyright (C) 2016-2022 The ESPResSo project
#
# This file is part of ESPResSo.
#
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

function(PYTHON_TEST)
  cmake_parse_arguments(TEST "" "FILE;MAX_NUM_PROC;SUFFIX"
                        "DEPENDS;DEPENDENCIES;LABELS;ARGUMENTS" ${ARGN})
  get_filename_component(TEST_NAME ${TEST_FILE} NAME_WE)
  set(TEST_FILE_CONFIGURED "${CMAKE_CURRENT_BINARY_DIR}/${TEST_NAME}.py")
  if(TEST_SUFFIX)
    set(TEST_NAME "${TEST_NAME}_${TEST_SUFFIX}")
  endif()
  configure_file(${TEST_FILE} ${TEST_FILE_CONFIGURED} COPYONLY)
  foreach(dependency IN LISTS TEST_DEPENDENCIES)
    configure_file(${dependency} ${CMAKE_CURRENT_BINARY_DIR}/${dependency}
                   COPYONLY)
  endforeach(dependency)
  set(TEST_FILE ${TEST_FILE_CONFIGURED})

  if(NOT DEFINED TEST_MAX_NUM_PROC)
    set(TEST_MAX_NUM_PROC 4)
  endif()

  if(${TEST_MAX_NUM_PROC} GREATER ${TEST_NP})
    set(TEST_NUM_PROC ${TEST_NP})
  else()
    set(TEST_NUM_PROC ${TEST_MAX_NUM_PROC})
  endif()

  if(EXISTS ${MPIEXEC})
    set_mpiexec_tmpdir("${TEST_NAME}")
    add_test(
      NAME ${TEST_NAME}
      COMMAND
        ${MPIEXEC} ${MPIEXEC_OVERSUBSCRIBE} ${MPIEXEC_NUMPROC_FLAG}
        ${TEST_NUM_PROC} ${MPIEXEC_PREFLAGS} ${MPIEXEC_TMPDIR}
        ${CMAKE_BINARY_DIR}/pypresso ${PYPRESSO_OPTIONS}
        ${TEST_FILE_CONFIGURED} ${TEST_ARGUMENTS} ${MPIEXEC_POSTFLAGS})
  else()
    add_test(${TEST_NAME} ${CMAKE_BINARY_DIR}/pypresso ${PYPRESSO_OPTIONS}
             ${TEST_FILE_CONFIGURED} ${TEST_ARGUMENTS})
  endif()
  set_tests_properties(${TEST_NAME} PROPERTIES PROCESSORS ${TEST_NUM_PROC}
                                               DEPENDS "${TEST_DEPENDS}")

  if("gpu" IN_LIST TEST_LABELS AND WITH_CUDA)
    set_tests_properties(${TEST_NAME} PROPERTIES RESOURCE_LOCK GPU)
  endif()

  if(${TEST_NUM_PROC} EQUAL 2)
    list(APPEND TEST_LABELS "parallel")
  endif()
  if(${TEST_NUM_PROC} EQUAL 3)
    list(APPEND TEST_LABELS "parallel_odd")
  endif()
  set_tests_properties(${TEST_NAME} PROPERTIES LABELS "${TEST_LABELS}")

  set(python_tests ${python_tests} ${TEST_FILE_CONFIGURED} PARENT_SCOPE)
endfunction(PYTHON_TEST)

function(CHECKPOINT_TEST)
  cmake_parse_arguments(TEST "" "MODES;MAX_NUM_PROC;SUFFIX" "LABELS" ${ARGN})
  if(NOT DEFINED TEST_MAX_NUM_PROC)
    set(TEST_MAX_NUM_PROC 4)
  endif()
  if(TEST_SUFFIX)
    set(TEST_ARGUMENTS "Test_suffix_${TEST_SUFFIX}__${TEST_MODES}")
    set(TEST_SUFFIX "_${TEST_MODES}_${TEST_SUFFIX}")
  else()
    set(TEST_ARGUMENTS "Test__${TEST_MODES}")
    set(TEST_SUFFIX "_${TEST_MODES}")
  endif()
  python_test(
    FILE save_checkpoint.py MAX_NUM_PROC ${TEST_MAX_NUM_PROC} LABELS
    ${TEST_LABELS} SUFFIX ${TEST_SUFFIX} ARGUMENTS ${TEST_ARGUMENTS}
    DEPENDENCIES unittest_generator.py)
  python_test(
    FILE
    test_checkpoint.py
    MAX_NUM_PROC
    ${TEST_MAX_NUM_PROC}
    LABELS
    ${TEST_LABELS}
    SUFFIX
    ${TEST_SUFFIX}
    ARGUMENTS
    ${TEST_ARGUMENTS}
    DEPENDENCIES
    unittest_generator.py
    DEPENDS
    save_checkpoint_${TEST_SUFFIX})
endfunction(CHECKPOINT_TEST)

python_test(FILE canonical_ensemble.py MAX_NUM_PROC 2)

add_custom_target(
  python_test_data
  COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/data
          ${CMAKE_CURRENT_BINARY_DIR}/data
  COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/tests_common.py
          ${CMAKE_CURRENT_BINARY_DIR}
  COMMAND
    ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/unittest_decorators.py
    ${CMAKE_CURRENT_BINARY_DIR}
  COMMAND
    ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/thermostats_common.py
    ${CMAKE_CURRENT_BINARY_DIR})

add_custom_target(
  check_python_parallel_odd
  COMMAND ${CMAKE_CTEST_COMMAND} --timeout ${TEST_TIMEOUT} -L parallel_odd
          ${CTEST_ARGS} --output-on-failure)
add_dependencies(check_python_parallel_odd pypresso python_test_data)
add_custom_target(
  check_python_gpu COMMAND ${CMAKE_CTEST_COMMAND} --timeout ${TEST_TIMEOUT} -L
                           gpu ${CTEST_ARGS} --output-on-failure)
add_dependencies(check_python_gpu pypresso python_test_data)

add_custom_target(
  check_python_skip_long
  COMMAND ${CMAKE_CTEST_COMMAND} --timeout ${TEST_TIMEOUT} -LE long
          ${CTEST_ARGS} --output-on-failure)
add_dependencies(check_python_skip_long pypresso python_test_data)

add_custom_target(
  check_python COMMAND ${CMAKE_CTEST_COMMAND} --timeout ${TEST_TIMEOUT}
                       ${CTEST_ARGS} --output-on-failure)
add_dependencies(check_python pypresso python_test_data)
add_dependencies(check check_python)
