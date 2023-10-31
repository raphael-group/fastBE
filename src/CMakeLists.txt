configure_file(
    "${fastbe_SOURCE_DIR}/include/fastbe.hpp.in"
    "${fastbe_SOURCE_DIR}/include/fastbe.hpp"
)

add_executable(fastbe fastbe.cxx)

# add libraries
target_link_libraries(fastbe PRIVATE nlohmann_json::nlohmann_json)
target_link_libraries(fastbe PRIVATE pprint)
target_link_libraries(fastbe PRIVATE spdlog)
target_link_libraries(fastbe PRIVATE csv)
target_link_libraries(fastbe PRIVATE argparse)
