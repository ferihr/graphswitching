/* Copyright (C) 2020--2023, Ferdinand Ihringer
 *
 * This file is part of Graph Switching.
 *
 * Graph Switching is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 */

#ifndef GRAPHSWITCHING_H
#define GRAPHSWITCHING_H

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

#define GRAPHSWITCHING_MAX_VERTICES 448
#define GRAPHSWITCHING_MAX_PART_SIZE 8

enum graphswitching_result {
        GRAPHSWITCHING_SUCCESS = 0,
        GRAPHSWITCHING_INVALID_ARGUMENT,
        GRAPHSWITCHING_INPUT_ERROR,
        GRAPHSWITCHING_OUTPUT_ERROR
};

/*
 * Read a vertex_count by vertex_count adjacency matrix from input and write
 * every graph obtainable by one WQH switching operation to output.
 *
 * The two switching parts both have part_size vertices. Input characters
 * other than '0' and '1' are ignored. The output starts with "n=<vertices>"
 * and uses adjacency-matrix format.
 */
enum graphswitching_result graphswitching_generate(
        FILE *input,
        FILE *output,
        int vertex_count,
        int part_size);

const char *graphswitching_result_string(enum graphswitching_result result);

#ifdef __cplusplus
}
#endif

#endif
