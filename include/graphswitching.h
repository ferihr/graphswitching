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

#define GRAPHSWITCHING_MAX_VERTICES 967
#define GRAPHSWITCHING_MAX_PART_SIZE 8
#define GRAPHSWITCHING_VERSION_MAJOR 0
#define GRAPHSWITCHING_VERSION_MINOR 1
#define GRAPHSWITCHING_VERSION_PATCH 7
#ifndef GRAPHSWITCHING_VERSION
#define GRAPHSWITCHING_VERSION "0.1.7"
#endif

enum graphswitching_method {
        GRAPHSWITCHING_METHOD_GM = 0,
        GRAPHSWITCHING_METHOD_WQH,
        GRAPHSWITCHING_METHOD_GM6,
        GRAPHSWITCHING_METHOD_WQH6,
        GRAPHSWITCHING_METHOD_AH6,
        GRAPHSWITCHING_METHOD_IS6,
        GRAPHSWITCHING_METHOD_FANO,
        GRAPHSWITCHING_METHOD_GM8,
        GRAPHSWITCHING_METHOD_GM44,
        GRAPHSWITCHING_METHOD_WQH8,
        GRAPHSWITCHING_METHOD_IS8_LEVEL3,
        GRAPHSWITCHING_METHOD_IS8_LEVEL5,
        GRAPHSWITCHING_METHOD_AH10
};

enum graphswitching_symmetry_mode {
        GRAPHSWITCHING_SYMMETRY_OFF = 0,
        GRAPHSWITCHING_SYMMETRY_ON,
        GRAPHSWITCHING_SYMMETRY_AUTO
};

enum graphswitching_output_format {
        GRAPHSWITCHING_OUTPUT_MATRIX = 0,
        GRAPHSWITCHING_OUTPUT_GRAPH6
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
        /*
         * Select ordinary, forced symmetry, or automatic symmetry search.
         * Both non-off modes require GRAPHSWITCHING_WITH_NAUTY.
         */
        enum graphswitching_symmetry_mode symmetry_mode;
        /* Deprecated compatibility alias: nonzero means forced symmetry. */
        int use_symmetry;
        /* Select adjacency-matrix or one-record-per-line graph6 output. */
        enum graphswitching_output_format output_format;
};

enum graphswitching_result {
        GRAPHSWITCHING_SUCCESS = 0,
        GRAPHSWITCHING_INVALID_ARGUMENT,
        GRAPHSWITCHING_INPUT_ERROR,
        GRAPHSWITCHING_OUTPUT_ERROR,
        GRAPHSWITCHING_FEATURE_UNAVAILABLE,
        GRAPHSWITCHING_MEMORY_ERROR
};

/* Initialize options for GM switching with automatic order detection. */
void graphswitching_options_init(struct graphswitching_options *options);

/*
 * Read one adjacency matrix and write every graph obtainable by one
 * switching operation. Matrix output starts with "n=<vertices>"; graph6
 * output contains one graph per line without a header.
 *
 * GM switching uses a four-vertex switching set. WQH switching uses two
 * parts of part_size vertices each. The remaining method values select the
 * fixed irreducible Simoens--Van Overberghe catalogues. With symmetry search
 * enabled, output multiplicity and order can change, but every result omitted
 * is isomorphic to a result produced from an orbit representative. Automatic
 * mode uses symmetry only when the computed automorphism group is large
 * enough to repay the orbit-search overhead.
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
