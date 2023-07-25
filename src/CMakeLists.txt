configure_file(
    "${vafpp_SOURCE_DIR}/include/vafpp.hpp.in"
    "${vafpp_SOURCE_DIR}/include/vafpp.hpp"
)

add_executable(vafpp vafpp.cxx)

# add libraries
target_link_libraries(vafpp PRIVATE nlohmann_json::nlohmann_json)
target_link_libraries(vafpp PRIVATE pprint)
target_link_libraries(vafpp PRIVATE spdlog)
target_link_libraries(vafpp PRIVATE csv)
target_link_libraries(vafpp PRIVATE argparse)
