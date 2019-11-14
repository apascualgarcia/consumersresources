cmake_minimum_required(VERSION 3.0)
project(consumerresources)

set(CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_INSTALL_LIBDIR additional_libs)
add_compile_options(-Wall -pedantic -O)

file(GLOB CRM_HEADERS include/*.h)
file(GLOB CRM_SOURCES src/CRModel_source/*.cc)
file(GLOB APP_SOURCES src/CRModel_targets/*.cc)

## we take care of installing libraries
# if we can’t find Eigen3, we install it
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
if(NOT Eigen3_FOUND)
	message(STATUS “Eigen3 not found, we will install it now.”)
	include(GNUInstallDirs)
	install(
		TARGETS
			Eigen3
		LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
		ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
		RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
	)
endif(NOT Eigen3_FOUND)

# if we can’t find gsl, we install it
find_package(GSL REQUIRED)
if(NOT GSL_FOUND)
	message(STATUS "GSL not found, we will install it now.")
	include(GNUInstallDirs)
	install(
		TARGETS
			GSL
		LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
		ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
		RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
	)
endif(NOT GSL_FOUND)

add_library(CRModel SHARED ${CRM_SOURCES} ${CRM_HEADERS})
target_include_directories(CRModel PUBLIC include)
target_link_libraries(CRModel Eigen3::Eigen GSL::gsl)

foreach(app ${APP_SOURCES})
	get_filename_component(app_name ${app} NAME_WE)
	add_executable(${app_name} ${app})
	target_link_libraries(${app_name} CRModel)
endforeach(app ${APP_SOURCES})
