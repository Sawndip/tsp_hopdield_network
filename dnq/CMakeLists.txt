cmake_minimum_required(VERSION 2.8) # Проверка версии CMake.
									# Если версия установленой программы
									# старее указаной, произайдёт аварийный выход.
set(CMAKE_${lang}_FLAGS_DEBUG_INIT "/Ox ${_RTC1}")

project(dnq)				# Название проекта

set(SOURCE_LIB dnq.cpp)		# Установка переменной со списком исходников
include_directories(nwta)			# Расположение заголовочных файлов

add_library(dnq STATIC ${SOURCE_LIB})# Создание статической библиотеки
find_package(Armadillo CONFIG REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})

target_link_libraries(dnq ${ARMADILLO_LIBRARIES})

target_link_libraries(dnq nwta)	# Линковка программы с библиотекой