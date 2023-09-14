#!/bin/bash
echo "Google Test started"
# Define your compile command
COMPILE_COMMAND="clang++ -o tests tests.cpp ../src/ErrorCorrection/hammingCode.cpp -I /usr/local/Cellar/googletest/1.14.0/include -lgtest -lgtest_main -pthread -std=c++14"
$COMPILE_COMMAND
# Define your test command
TEST_COMMAND="./tests"
$TEST_COMMAND