project(multiStateFlowReliability)
cmake_minimum_required(VERSION 3.1)
set(Boost_USE_STATIC_LIBS	OFF)
set(Boost_USE_MULTITHREADED	ON)
set(Boost_USE_STATIC_RUNTIME	OFF)

if(CMAKE_COMPILER_IS_GNUCC)
	add_definitions("-Wall")
endif()

#Find boost packages
find_package(Boost 1.53.0 COMPONENTS graph program_options serialization regex REQUIRED)

if(WIN32)
	find_package(mpfr)
	find_package(mpir)

	#If we're on windows, Start copying dlls
	set(BUILD_TYPES Release Debug RelWithDebInfo MinSizeRel)
	set(BOOST_BUILD_TYPES Release Debug Release Release)
	foreach(val RANGE 3)
		list(GET BUILD_TYPES ${val} BUILD_TYPE)
		list(GET BOOST_BUILD_TYPES ${val} BOOST_BUILD_TYPE)
		string(TOUPPER ${BOOST_BUILD_TYPE} BOOST_BUILD_TYPE_UPPER)
		
		#copy boost serialization
		get_filename_component(NAME_PART ${Boost_SERIALIZATION_LIBRARY_${BOOST_BUILD_TYPE_UPPER}} NAME_WE)
		get_filename_component(PATH_PART ${Boost_SERIALIZATION_LIBRARY_${BOOST_BUILD_TYPE_UPPER}} PATH)
		file(COPY ${PATH_PART}/${NAME_PART}.dll DESTINATION ${CMAKE_BINARY_DIR}/${BUILD_TYPE}/)

		#copy boost program_options
		get_filename_component(NAME_PART ${Boost_PROGRAM_OPTIONS_LIBRARY_${BOOST_BUILD_TYPE_UPPER}} NAME_WE)
		get_filename_component(PATH_PART ${Boost_PROGRAM_OPTIONS_LIBRARY_${BOOST_BUILD_TYPE_UPPER}} PATH)
		file(COPY ${PATH_PART}/${NAME_PART}.dll DESTINATION ${CMAKE_BINARY_DIR}/${BUILD_TYPE}/)

		#copy boost graph
		get_filename_component(NAME_PART ${Boost_GRAPH_LIBRARY_${BOOST_BUILD_TYPE_UPPER}} NAME_WE)
		get_filename_component(PATH_PART ${Boost_GRAPH_LIBRARY_${BOOST_BUILD_TYPE_UPPER}} PATH)
		file(COPY ${PATH_PART}/${NAME_PART}.dll DESTINATION ${CMAKE_BINARY_DIR}/${BUILD_TYPE}/)

		#copy boost regex
		get_filename_component(NAME_PART ${Boost_REGEX_LIBRARY_${BOOST_BUILD_TYPE_UPPER}} NAME_WE)
		get_filename_component(PATH_PART ${Boost_REGEX_LIBRARY_${BOOST_BUILD_TYPE_UPPER}} PATH)
		file(COPY ${PATH_PART}/${NAME_PART}.dll DESTINATION ${CMAKE_BINARY_DIR}/${BUILD_TYPE}/)

		#Now for MPFR
		get_target_property(MPFR_LOCATION mpfr "LOCATION_${BOOST_BUILD_TYPE}")
		file(COPY "${MPFR_LOCATION}" DESTINATION "${CMAKE_BINARY_DIR}/${BUILD_TYPE}/")

		#Now for MPIR
		get_target_property(MPIR_LOCATION mpir "LOCATION_${BOOST_BUILD_TYPE}")
		file(COPY "${MPIR_LOCATION}" DESTINATION "${CMAKE_BINARY_DIR}/${BUILD_TYPE}/")
	endforeach()
else()
	#mpfr needs to be listed as an imported target so that we can use target_link_librarie later on. If we don't have this then it just adds a -lmpfr which probably picks up the wrong thing. 
	find_library(mpfr_LIBRARY mpfr)
	add_library(mpfr SHARED IMPORTED)
	set_target_properties(mpfr PROPERTIES IMPORTED_LOCATION ${mpfr_LIBRARY})
endif()
add_subdirectory(common)
add_subdirectory(crudeMC)
add_subdirectory(PMCExecutable)
add_subdirectory(turnip)
