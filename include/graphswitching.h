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
#define GRAPHSWITCHING_VERSION_MAJOR 0
#define GRAPHSWITCHING_VERSION_MINOR 1
#define GRAPHSWITCHING_VERSION_PATCH 1
#ifndef GRAPHSWITCHING_VERSION
#define GRAPHSWITCHING_VERSION "0.1.1"
#endif

enum graphswitching_method {
        GRAPHSWITCHING_METHOD_GM = 0,
        GRAPHSWITCHING_METHOD_WQH
};

struct graphswitching_options {
        enum graphswitching_method method;
        /*
         * Set vertex_count to zero to infer the order from a square matrix
         * or from an optional leading "n=<vertices>" line.
         */
        int vertex_count;
        /* Used only by WQH switching. */
        int part_size;
};

enum graphswitching_result {
        GRAPHSWITCHING_SUCCESS = 0,
        GRAPHSWITCHING_INVALID_ARGUMENT,
        GRAPHSWITCHING_INPUT_ERROR,
        GRAPHSWITCHING_OUTPUT_ERROR
};

/* Initialize options for GM switching with automatic order detection. */
void graphswitching_options_init(struct graphswitching_options *options);

/*
 * Read one adjacency matrix and write every graph obtainable by one
 * switching operation. The output starts with "n=<vertices>" and uses
 * adjacency-matrix format.
 *
 * GM switching uses a four-vertex switching set. WQH switching uses two
 * parts of part_size vertices each.
 */
enum graphswitching_result graphswitching_generate_with_options(
        FILE *input,
        FILE *output,
        const struct graphswitching_options *options);

/*
 * Compatibility interface for WQH switching. New callers should use
 * graphswitching_generate_with_options().
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
