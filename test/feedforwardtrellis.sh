#!/bin/bash
echo "Google Test started"
# Define your compile command
COMPILE_COMMAND="clang++ -g -o trellis_test feedforwardtrellis_test.cpp ../src/ErrorCorrection/feedForwardTrellis.cpp -I /usr/local/Cellar/googletest/1.14.0/include -I ../src/helper.h -lgtest -lgtest_main -pthread -std=c++14"
$COMPILE_COMMAND
# Define your test command
TEST_COMMAND="./trellis_test"
$TEST_COMMAND