# Copyright (c) 2023 Jim-shop
# bp is licensed under Mulan PubL v2.
# You can use this software according to the terms and conditions of the Mulan PubL v2.
# You may obtain a copy of Mulan PubL v2 at:
#          http://license.coscl.org.cn/MulanPubL-2.0
# THIS SOFTWARE IS PROVIDED ON AN "AS IS" BASIS, WITHOUT WARRANTIES OF ANY KIND,
# EITHER EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO NON-INFRINGEMENT,
# MERCHANTABILITY OR FIT FOR A PARTICULAR PURPOSE.
# See the Mulan PubL v2 for more details.
cmake_minimum_required(VERSION 3.12)

project(bp
    VERSION 0.1.0.0
    DESCRIPTION "belief propagation decoder"
    HOMEPAGE_URL https://github.com/Jim-shop/bp
    LANGUAGES CXX
)

set(CMAKE_BUILD_TYPE Release)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_C_STANDARD 17)

add_library(SparseMatrix STATIC
    src/lib/sparse_matrix/sparse_matrix.cpp
)
target_include_directories(SparseMatrix
    PUBLIC src/lib/sparse_matrix
)

add_library(Json INTERFACE)
target_include_directories(Json
    INTERFACE src/lib/nlohmann/
)

add_library(BpDecoder STATIC
    src/lib/bp_decoder/bp_decoder.cpp
)
target_include_directories(BpDecoder
    PUBLIC src/lib/bp_decoder
    PUBLIC src/lib/sparse_matrix
)
target_link_libraries(BpDecoder
    PRIVATE SparseMatrix
)

add_executable(sim
    src/sim/sim.cpp
)
target_link_libraries(sim
    PRIVATE SparseMatrix
    PRIVATE Json
    PRIVATE BpDecoder
)
target_include_directories(sim
    PRIVATE src/lib/
)