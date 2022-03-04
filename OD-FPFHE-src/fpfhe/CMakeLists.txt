#
# CMakeLists.txt for FPFHE library
#

# all files named *.cpp are compiled to form the library
file (GLOB FPFHE_SRC_FILES CONFIGURE_DEPENDS lib/*.cpp)

include_directories(${CORE_INCLUDE_DIRS})
list(APPEND FPFHE_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include")
list(APPEND FPFHE_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/lib")
include_directories(${FPFHE_INCLUDE_DIRS})

set(FPFHE_VERSION_MAJOR ${PALISADE_VERSION_MAJOR})
set(FPFHE_VERSION_MINOR ${PALISADE_VERSION_MINOR})
set(FPFHE_VERSION_PATCH ${PALISADE_VERSION_PATCH})
set(FPFHE_VERSION ${FPFHE_VERSION_MAJOR}.${FPFHE_VERSION_MINOR}.${FPFHE_VERSION_PATCH})

add_library(fpfheobj OBJECT ${FPFHE_SRC_FILES})
set_property(TARGET fpfheobj PROPERTY POSITION_INDEPENDENT_CODE 1)

if ( BUILD_SHARED )
	add_dependencies(fpfheobj PALISADEcore)
	add_library (PALISADEfpfhe SHARED $<TARGET_OBJECTS:fpfheobj>)
	set_property(TARGET PALISADEfpfhe PROPERTY VERSION ${FPFHE_VERSION})
	set_property(TARGET PALISADEfpfhe PROPERTY SOVERSION ${FPFHE_VERSION_MAJOR})
	set_property(TARGET PALISADEfpfhe PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
	install(TARGETS PALISADEfpfhe
		EXPORT PalisadeTargets
		DESTINATION lib)
endif()

if( BUILD_STATIC )
add_dependencies(fpfheobj PALISADEcore_static)
	add_library (PALISADEfpfhe_static STATIC $<TARGET_OBJECTS:fpfheobj>)
	set_property(TARGET PALISADEfpfhe_static PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
	install(TARGETS PALISADEfpfhe_static
		EXPORT PalisadeTargets
		DESTINATION lib)
endif()

install(DIRECTORY include/	DESTINATION include/palisade/fpfhe)

add_custom_target( allfpfhe )

if( BUILD_SHARED )
	set (FPFHELIBS PUBLIC PALISADEfpfhe PUBLIC PALISADEcore ${THIRDPARTYLIBS} ${OpenMP_CXX_FLAGS})
	target_link_libraries (PALISADEfpfhe PUBLIC PALISADEcore ${THIRDPARTYLIBS} ${OpenMP_CXX_FLAGS})
	add_dependencies( allfpfhe PALISADEfpfhe  )
endif()

if( BUILD_STATIC )
	set (FPFHELIBS ${FPFHELIBS} PUBLIC PALISADEfpfhe_static PUBLIC PALISADEcore_static ${THIRDPARTYSTATICLIBS} ${OpenMP_CXX_FLAGS})
	target_link_libraries (PALISADEfpfhe_static PUBLIC PALISADEcore_static ${THIRDPARTYSTATICLIBS} ${OpenMP_CXX_FLAGS})
	add_dependencies( allfpfhe PALISADEfpfhe_static )
endif()

if( BUILD_UNITTESTS )
	# UNITEST Later~! SH

#	file (GLOB FPFHE_TEST_SRC_FILES CONFIGURE_DEPENDS unittest/*.cpp)
#	add_executable (fpfhe_tests ${BINFHE_TEST_SRC_FILES} ${UNITTESTMAIN})
#	set_property(TARGET fpfhe_tests PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/unittest)
#	target_link_libraries ( fpfhe_tests ${BINFHELIBS})
#	if (NOT ${WITH_OPENMP})
#	    target_link_libraries ( fpfhe_tests PRIVATE Threads::Threads)
#	endif()

#	add_dependencies( allfpfhe fpfhe_tests )

#	add_custom_command( OUTPUT runfpfhetests WORKING_DIRECTORY ${CMAKE_BINARY_DIR} COMMAND
#   ${CMAKE_BINARY_DIR}/unittest/fpfhe_tests )
#	add_custom_target( testfpfhe DEPENDS fpfhe_tests runfpfhetests )
endif()




# Example later!! SH
#set (FPFHEAPPS "")
#if( BUILD_EXAMPLES)
#    file (GLOB FPFHE_EXAMPLES_SRC_FILES CONFIGURE_DEPENDS examples/*.cpp)
#	foreach (app ${FPFHE_EXAMPLES_SRC_FILES})
#        get_filename_component ( exe ${app} NAME_WE )
#		add_executable ( ${exe} ${app} ) 
#		set_property(TARGET ${exe} PROPERTY RUNTIME_OUTPUT_DIRECTORY
#        ${CMAKE_BINARY_DIR}/bin/examples/fpfhe)
#		set( FPFHEAPPS ${FPFHEAPPS} ${exe} )
#		target_link_libraries ( ${exe} ${FPFHELIBS} )
#	endforeach()

#	add_custom_target( allfpfheexamples )
	#add_dependencies( allfpfheexamples ${FPFHEAPPS} )
	#add_dependencies( allfpfhe allfpfheexamples)
#endif()

#add_custom_command( OUTPUT fpfheinfocmd COMMAND echo Builds PALISADEfpfhe and these apps:${FPFHEAPPS} )
#add_custom_target( fpfheinfo DEPENDS fpfheinfocmd )
