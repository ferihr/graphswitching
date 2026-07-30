/* Copyright (C) 2020--2023, Ferdinand Ihringer
 *
 * Command-line interface for generic WQH graph switching.
 */

#include "graphswitching.h"

#include <errno.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>

#define OUTPUT_BUFFER_SIZE (1024U * 1024U)

static int parse_positive_int(const char *text, int *value)
{
        char *end;
        long parsed;

        errno = 0;
        end = NULL;
        parsed = strtol(text, &end, 10);
        if (errno == ERANGE || end == text || *end != '\0' ||
            parsed <= 0 || parsed > INT_MAX) {
                return 0;
        }

        *value = (int)parsed;
        return 1;
}

static void print_usage(FILE *stream, const char *program)
{
        fprintf(stream, "Usage: %s n p\n", program);
        fprintf(stream,
                "Read an n-by-n adjacency matrix and apply WQH switching "
                "with parts of size p.\n");
}

int main(int argc, char *argv[])
{
        static char output_buffer[OUTPUT_BUFFER_SIZE];
        enum graphswitching_result result;
        int vertex_count;
        int part_size;

        /*
         * Matrix enumeration can produce a large amount of output. A larger
         * full buffer reduces the number of writes without changing the
         * output format.
         */
        (void)setvbuf(stdout, output_buffer, _IOFBF, sizeof(output_buffer));

        if (argc != 3) {
                print_usage(stderr, argv[0]);
                return EXIT_FAILURE;
        }

        if (!parse_positive_int(argv[1], &vertex_count) ||
            !parse_positive_int(argv[2], &part_size)) {
                fprintf(stderr, "%s: n and p must be positive integers\n",
                        argv[0]);
                return EXIT_FAILURE;
        }

        if (vertex_count > GRAPHSWITCHING_MAX_VERTICES) {
                fprintf(stderr, "%s: n must not exceed %d\n", argv[0],
                        GRAPHSWITCHING_MAX_VERTICES);
                return EXIT_FAILURE;
        }

        if (part_size > GRAPHSWITCHING_MAX_PART_SIZE) {
                fprintf(stderr, "%s: p must not exceed %d\n", argv[0],
                        GRAPHSWITCHING_MAX_PART_SIZE);
                return EXIT_FAILURE;
        }

        if (2 * part_size > vertex_count) {
                fprintf(stderr, "%s: n must be at least 2p\n", argv[0]);
                return EXIT_FAILURE;
        }

        result = graphswitching_generate(stdin, stdout, vertex_count,
                                         part_size);
        if (result != GRAPHSWITCHING_SUCCESS) {
                fprintf(stderr, "%s: %s\n", argv[0],
                        graphswitching_result_string(result));
                return EXIT_FAILURE;
        }

        if (fflush(stdout) == EOF) {
                fprintf(stderr, "%s: could not write output\n", argv[0]);
                return EXIT_FAILURE;
        }

        return EXIT_SUCCESS;
}
