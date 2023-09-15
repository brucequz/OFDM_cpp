#!/bin/bash
echo "Google Test started"
# Define your compile command
COMPILE_COMMAND="clang++ -o hamming_test hamming_test.cpp ../src/ErrorCorrection/hammingCode.cpp -I /usr/local/Cellar/googletest/1.14.0/include -lgtest -lgtest_main -pthread -std=c++14"
$COMPILE_COMMAND
# Define your test command
TEST_COMMAND="./hamming_test"
$TEST_COMMAND