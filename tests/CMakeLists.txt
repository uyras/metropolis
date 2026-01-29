
add_test(
  NAME
    ising
  COMMAND
    $<TARGET_FILE:metropolis> ${CMAKE_SOURCE_DIR}/examples/ising/ising_16_simple.ini
)
