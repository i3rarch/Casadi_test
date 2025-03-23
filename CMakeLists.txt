cmake_minimum_required(VERSION 3.10)
project(CasadiExample)

set(CMAKE_CXX_STANDARD 14)

# Путь к директории CasADi 
set(CASADI_DIR $ENV{USERPROFILE}/AppData/Local/Programs/Python/Python313/Lib/site-packages/casadi)

# Добавляем исполняемый файл
add_executable(casadi_example main.cpp)

# Подключаем заголовочные файлы CasADi
target_include_directories(casadi_example PRIVATE ${CASADI_DIR}/include)

# Настройки для Windows
if(WIN32)
    # Ищем библиотеку CasADi
    find_library(CASADI_LIB 
        NAMES casadi libcasadi
        PATHS ${CASADI_DIR}
        NO_DEFAULT_PATH
    )
    
    if(NOT CASADI_LIB)
        message(FATAL_ERROR "Библиотека CasADi не найдена в ${CASADI_DIR}")
    endif()
    
    message(STATUS "Найдена библиотека CasADi: ${CASADI_LIB}")
    
    # Линкуем только основную библиотеку
    target_link_libraries(casadi_example ${CASADI_LIB})
    
    # Добавляем путь к DLL для запуска
    set_target_properties(casadi_example PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
        VS_DEBUGGER_ENVIRONMENT "PATH=${CASADI_DIR};$ENV{PATH}"
    )

    # Копируем DLL в директорию сборки
    file(GLOB CASADI_DLLS "${CASADI_DIR}/*.dll")
    foreach(DLL ${CASADI_DLLS})
        add_custom_command(TARGET casadi_example POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
            ${DLL}
            ${CMAKE_BINARY_DIR})
    endforeach()
endif()