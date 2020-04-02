cmake_minimum_required(VERSION 2.8) # Проверка версии CMake.
									# Если версия установленой программы
									# старее указаной, произайдёт аварийный выход.
set(CMAKE_${lang}_FLAGS_DEBUG_INIT "/Ox ${_RTC1}")

project(nwta)				# Название проекта

set(SOURCE_LIB nwta.cpp 2opt/2opt.h 2opt/2opt.cpp)		# Установка переменной со списком исходников

add_library(nwta STATIC ${SOURCE_LIB})# Создание статической библиотеки