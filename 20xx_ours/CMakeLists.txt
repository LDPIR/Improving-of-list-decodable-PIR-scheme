cmake_minimum_required(VERSION 3.0.0)

project(20xx_ours VERSION 0.1.0)

# 添加源文件和头文件
file(GLOB SOURCES src/*.cpp)
file(GLOB HEADERS include/*.hpp)

# 输出源文件和头文件列表（可选）
# message("Source files: ${SOURCES}")
# message("Header files: ${HEADERS}")

# 添加可执行文件
add_executable(20xx_ours ${SOURCES} ${HEADERS})

# 查找并链接 libgmp 库
find_library(GMP_LIB gmp)
if (GMP_LIB)
    target_link_libraries(20xx_ours PRIVATE ${GMP_LIB})
else()
    message(FATAL_ERROR "libgmp not found")
endif()

target_link_libraries(20xx_ours PRIVATE m gmp flint)


# 设置包含路径（如果有必要）
target_include_directories(20xx_ours PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

# 设置编译定义（如果有必要）
target_compile_definitions(20xx_ours PRIVATE DATA_DIR="${DATA_DIR}")
