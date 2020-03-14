cmake_minimum_required(VERSION 2.8) # Проверка версии CMake.
									# Если версия установленой программы
									# старее указаной, произайдёт аварийный выход.

project(nwta)				# Название проекта

set(SOURCE_LIB nwta.cpp)		# Установка переменной со списком исходников

add_library(nwta STATIC ${SOURCE_LIB})# Создание статической библиотеки