cmake_minimum_required(VERSION 2.8) # Проверка версии CMake.
									# Если версия установленой программы
									# старее указаной, произайдёт аварийный выход.
set(CMAKE_${lang}_FLAGS_DEBUG_INIT "/Ox ${_RTC1}")

project(hello_world)				# Название проекта

set(SOURCE_EXE main.cpp)			# Установка переменной со списком исходников

include_directories(nwta)			# Расположение заголовочных файлов
include_directories(dnq)			# Расположение заголовочных файлов

add_executable(main ${SOURCE_EXE})	# Создает исполняемый файл с именем main

add_subdirectory(nwta)				# Добавление подпроекта, указывается имя дирректории
add_subdirectory(dnq)				# Добавление подпроекта, указывается имя дирректории

find_package(Armadillo CONFIG REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})


target_link_libraries(main nwta dnq)	# Линковка программы с библиотекой
target_link_libraries(main ${ARMADILLO_LIBRARIES})
