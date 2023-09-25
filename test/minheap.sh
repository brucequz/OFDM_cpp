#!/bin/bash
echo "Google Test started"
# Define your remove command
REMOVE_COMMAND="rm minheap_test"
$REMOVE_COMMAND
# Define your compile command
COMPILE_COMMAND="clang++ -g -o minheap_test minheap_test.cpp ../src/ErrorCorrection/minHeap.cpp -I /usr/local/Cellar/googletest/1.14.0/include -lgtest -lgtest_main -pthread -std=c++14"
$COMPILE_COMMAND
# Define your test command
TEST_COMMAND="./minheap_test"
$TEST_COMMAND