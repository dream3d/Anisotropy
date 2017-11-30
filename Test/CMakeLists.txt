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


list(APPEND ${PLUGIN_NAME}_LINK_LIBS Qt5::Core H5Support SIMPLib)

if(APPLE)
 list(APPEND ${PLUGIN_NAME}_LINK_LIBS SimpleITKExplicit)
else()
  list(APPEND ${PLUGIN_NAME}_LINK_LIBS ${ITK_LIBRARIES})
endif()

AddSIMPLUnitTest(TESTNAME ${PLUGIN_NAME}UnitTest
                  SOURCES ${${PLUGIN_NAME}_SOURCE_DIR}/Test/AnisotropyFilterTest.cpp
                  FOLDER "${PLUGIN_NAME}Plugin/Test"
                  LINK_LIBRARIES ${${PLUGIN_NAME}_LINK_LIBS}
                  INCLUDE_DIRS
                          ${PLUGINS_SOURCE_DIR}
                          ${PLUGINS_BINARY_DIR}
                          ${PLUGINS_BINARY_DIR}/${PLUGIN_NAME}/Test
)

if(MSVC)
  set_source_files_properties(${${PLUGIN_NAME}Test_BINARY_DIR}/${PLUGIN_NAME}UnitTest.cpp PROPERTIES COMPILE_FLAGS /bigobj)
endif()
