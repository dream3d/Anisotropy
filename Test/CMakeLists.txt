# ============================================================================
# Copyright (c) 2016 Czech Academy of Sciences, Institute of Physics,
# Group of Bulk Nanomaterials and Interfaces, http://ams.fzu.cz
#
# Redistribution and use in source and binary forms, with or without modification,
# are permitted provided that the following conditions are met:
#
# Redistributions of source code must retain the above copyright notice, this
# list of conditions and the following disclaimer.
#
# Redistributions in binary form must reproduce the above copyright notice, this
# list of conditions and the following disclaimer in the documentation and/or
# other materials provided with the distribution.
#
# Neither the name of the Czech Academy of Sciences, nor the names of its 
# contributors may be used to endorse or promote products derived from this 
# software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
# USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# The code contained herein was partially funded by the followig grants:
#    Czech Science Foundation (GA CR), project no. GBP108/12/G043
#    Czech Ministry of Education, Youth and Sports (MSMT), project no. LM2015087
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

set(${PLUGIN_NAME}Test_SOURCE_DIR ${${PLUGIN_NAME}_SOURCE_DIR}/Test)
set(${PLUGIN_NAME}Test_BINARY_DIR ${${PLUGIN_NAME}_BINARY_DIR}/Test)

#------------------------------------------------------------------------------
# Include this file from the CMP Project
include(${CMP_SOURCE_DIR}/cmpCMakeMacros.cmake)
include(${SIMPLProj_SOURCE_DIR}/Source/SIMPLib/SIMPLibMacros.cmake)

set(TEST_TEMP_DIR ${${PLUGIN_NAME}Test_BINARY_DIR}/Temp)
# Make sure the directory is created during CMake time
file(MAKE_DIRECTORY ${TEST_TEMP_DIR})

set(TEST_DATA_DIR ${${PLUGIN_NAME}_SOURCE_DIR}/Data/${PLUGIN_NAME})
# Make sure the directory is created during CMake time
file(MAKE_DIRECTORY ${TEST_DATA_DIR})

file(GLOB templateFiles RELATIVE "${${PLUGIN_NAME}_SOURCE_DIR}/Data" "${${PLUGIN_NAME}_SOURCE_DIR}/Data/*")
    foreach(templateFile ${templateFiles})
	    set(srcTemplatePath ${${PLUGIN_NAME}_SOURCE_DIR}/Data/${templateFile})
        if(NOT IS_DIRECTORY ${srcTemplatePath})
           configure_file(
                    ${${PLUGIN_NAME}_SOURCE_DIR}/Data/${templateFile}
                    ${${PLUGIN_NAME}_BINARY_DIR}/Test/Data/${templateFile}
                    COPYONLY)
		endif(NOT IS_DIRECTORY ${srcTemplatePath})
    endforeach(templateFile)	

file(GLOB templateFiles RELATIVE "${${PLUGIN_NAME}_SOURCE_DIR}/Data/tif" "${${PLUGIN_NAME}_SOURCE_DIR}/Data/tif/*")
    foreach(templateFile ${templateFiles})
	    set(srcTemplatePath ${${PLUGIN_NAME}_SOURCE_DIR}/Data/tif/${templateFile})
        if(NOT IS_DIRECTORY ${srcTemplatePath})
           configure_file(
                    ${${PLUGIN_NAME}_SOURCE_DIR}/Data/tif/${templateFile}
                    ${${PLUGIN_NAME}_BINARY_DIR}/Test/Data/${templateFile}
                    COPYONLY)
		endif(NOT IS_DIRECTORY ${srcTemplatePath})
    endforeach(templateFile)	

configure_file(${${PLUGIN_NAME}Test_SOURCE_DIR}/TestFileLocations.h.in
               ${${PLUGIN_NAME}Test_BINARY_DIR}/${PLUGIN_NAME}TestFileLocations.h @ONLY IMMEDIATE)

configure_file(${SIMPLProj_SOURCE_DIR}/Resources/UnitTestSupport.hpp
               ${${PLUGIN_NAME}Test_BINARY_DIR}/UnitTestSupport.hpp COPYONLY IMMEDIATE)


list(APPEND ${PLUGIN_NAME}_LINK_LIBS Qt5::Core H5Support SIMPLib ${ITK_LIBRARIES})

AddSIMPLUnitTest(TESTNAME ${PLUGIN_NAME}UnitTest
                  SOURCES ${${PLUGIN_NAME}_SOURCE_DIR}/Test/AnisotropyFilterTest.cpp
                  FOLDER "${PLUGIN_NAME}Plugin/Test"
                  LINK_LIBRARIES ${${PLUGIN_NAME}_LINK_LIBS}
                  INCLUDE_DIRS
                          ${PLUGINS_SOURCE_DIR}
                          ${PLUGINS_BINARY_DIR}
                          ${PLUGINS_BINARY_DIR}/${PLUGIN_NAME}/Test
)

#------------------------------------------------------------------------------
# If Python is enabled, then enable the Python unit tests for this plugin
if(SIMPL_ENABLE_PYTHON)
    get_property(SIMPL_ANACONDA_OUTPUT_DIR GLOBAL PROPERTY SIMPL_ANACONDA_OUTPUT_DIR)
    get_property(SIMPL_PYTHON_TEST_SCRIPT GLOBAL PROPERTY SIMPL_PYTHON_TEST_SCRIPT)
    get_property(PYTHON_SITE_PACKAGES_NAME GLOBAL PROPERTY PYTHON_SITE_PACKAGES_NAME)

    #------------------------------
    # Initialize the PYTHON_TEST_INPUT_DIR variable to point to the "Testing Directory"
    file(TO_CMAKE_PATH "${${PLUGIN_NAME}Test_SOURCE_DIR}/Python" PYTHON_TEST_INPUT_DIR)
    #------------------------------
    # These names should match the names "EXACTLY" (including capitalization).
    # NO Spaces in the names (which means no spaces in the variable names)
    set(PLUGIN_PYTHON_TESTS
      # 01_Adaptive_Alignment_Import_Data # The conversion is already complete so this does not need to run
      02_Adaptive_Alignment_Misorientation_Zero_Shifts
      03_Adaptive_Alignment_Mutual_Information_SEM_Images
      04_Adaptive_Alignment_Steiner_Compact
    )

    set(TEST_SCRIPT_FILE_EXT "sh")
    if(WIN32)
      set(TEST_SCRIPT_FILE_EXT "bat")
    endif()
    foreach(test ${PLUGIN_PYTHON_TESTS})
        set(SIMPL_PYTHON_TEST_DRIVER "${SIMPL_ANACONDA_OUTPUT_DIR}/UnitTest/${PLUGIN_NAME}/${test}.${TEST_SCRIPT_FILE_EXT}")

        configure_file(${SIMPL_PYTHON_TEST_SCRIPT}
                        "${SIMPL_PYTHON_TEST_DRIVER}" @ONLY)

        add_test(NAME PY_${PLUGIN_NAME}_${test} COMMAND "${SIMPL_PYTHON_TEST_DRIVER}" )
    endforeach(test ${PLUGIN_PYTHON_TESTS})
    #------------------------------
    # Also setup a unit test for the base python unit test file that is generated as part
    # of the Pybind11 generated codes
    set(PYTHON_TEST_INPUT_DIR "${SIMPL_ANACONDA_OUTPUT_DIR}/UnitTest/${PLUGIN_NAME}")
    set(SIMPL_PYTHON_TEST_DRIVER "${SIMPL_ANACONDA_OUTPUT_DIR}/UnitTest/${PLUGIN_NAME}/${PLUGIN_NAME}.sh")
    set(test "${PLUGIN_NAME}_UnitTest")
    configure_file(${SIMPL_PYTHON_TEST_SCRIPT}  "${SIMPL_PYTHON_TEST_DRIVER}" @ONLY)
    add_test(NAME PY_${PLUGIN_NAME}_UnitTest COMMAND "${SIMPL_PYTHON_TEST_DRIVER}" )
endif()

