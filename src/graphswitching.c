/* Copyright (C) 2020--2023, Ferdinand Ihringer
 *
 * Generic Wang--Qiu--Hu graph switching.
 *
 * This implementation retains the algorithm and adjacency-matrix output of
 * code/gen_all_srgs_wqh_generic.c while providing a reusable interface.
 */

#include "graphswitching.h"

#include <stdint.h>
#include <string.h>

#define WORD_BITS 32
#define BLOCK_SHIFT 5
#define BLOCK_MASK (WORD_BITS - 1)
#define BLOCK_COUNT \
        ((GRAPHSWITCHING_MAX_VERTICES + WORD_BITS - 1) / WORD_BITS)

typedef uint32_t word_t;

struct partition_vertex {
        int vertex;
        int row_offset;
        int word_index;
        int bit_index;
};

struct search_context {
        int vertex_count;
        int part_size;
        word_t adjacency[GRAPHSWITCHING_MAX_VERTICES * BLOCK_COUNT];
        struct partition_vertex
                partition[2 * GRAPHSWITCHING_MAX_PART_SIZE];
        FILE *output;
};

static int partition_contains(const struct partition_vertex partition[],
                              int last, int vertex);
static int adjacency_bit(const struct search_context *context, int row,
                         int column);
static void apply_wqh(struct search_context *context,
                      const int neighbour_counts[]);
static enum graphswitching_result choose_partition(
        struct search_context *context, int current);
static enum graphswitching_result apply_partition_global(
        struct search_context *context);
static int test_regular_11(const struct search_context *context, int current);
static int test_regular_12(const struct search_context *context, int current);
static int test_regular_21(const struct search_context *context, int current);
static int test_regular_22(const struct search_context *context, int current);
static void get_degrees_11(const struct search_context *context, int current,
                           int degrees[]);
static void get_degrees_12(const struct search_context *context, int current,
                           int degrees[]);
static void get_degrees_21(const struct search_context *context, int current,
                           int degrees[]);
static void get_degrees_22(const struct search_context *context, int current,
                           int degrees[]);
static enum graphswitching_result read_adjacency_matrix(
        struct search_context *context, FILE *input);
static enum graphswitching_result write_adjacency_matrix(
        const struct search_context *context);

enum graphswitching_result graphswitching_generate(
        FILE *input,
        FILE *output,
        int vertex_count,
        int part_size)
{
        struct search_context context;
        enum graphswitching_result result;

        if (input == NULL || output == NULL || vertex_count <= 0 ||
            vertex_count > GRAPHSWITCHING_MAX_VERTICES || part_size <= 0 ||
            part_size > GRAPHSWITCHING_MAX_PART_SIZE ||
            2 * part_size > vertex_count) {
                return GRAPHSWITCHING_INVALID_ARGUMENT;
        }

        memset(&context, 0, sizeof(context));
        context.vertex_count = vertex_count;
        context.part_size = part_size;
        context.output = output;

        result = read_adjacency_matrix(&context, input);
        if (result != GRAPHSWITCHING_SUCCESS) {
                return result;
        }

        if (fprintf(output, "n=%d\n", vertex_count) < 0) {
                return GRAPHSWITCHING_OUTPUT_ERROR;
        }

        return choose_partition(&context, 0);
}

const char *graphswitching_result_string(enum graphswitching_result result)
{
        switch (result) {
        case GRAPHSWITCHING_SUCCESS:
                return "success";
        case GRAPHSWITCHING_INVALID_ARGUMENT:
                return "invalid arguments";
        case GRAPHSWITCHING_INPUT_ERROR:
                return "could not read a complete adjacency matrix";
        case GRAPHSWITCHING_OUTPUT_ERROR:
                return "could not write output";
        default:
                return "unknown error";
        }
}

static enum graphswitching_result read_adjacency_matrix(
        struct search_context *context, FILE *input)
{
        int row;
        int column;

        for (row = 0; row < context->vertex_count; ++row) {
                int row_offset = row * BLOCK_COUNT;

                for (column = 0; column < context->vertex_count; ++column) {
                        int character;

                        do {
                                character = fgetc(input);
                                if (character == EOF) {
                                        return GRAPHSWITCHING_INPUT_ERROR;
                                }
                        } while (character != '0' && character != '1');

                        if (character == '1') {
                                int word_index = column >> BLOCK_SHIFT;
                                int bit_index = column & BLOCK_MASK;

                                context->adjacency[row_offset + word_index] |=
                                        UINT32_C(1) << bit_index;
                        }
                }
        }

        return GRAPHSWITCHING_SUCCESS;
}

static int adjacency_bit(const struct search_context *context, int row,
                         int column)
{
        int word_index = column >> BLOCK_SHIFT;
        int bit_index = column & BLOCK_MASK;

        return (int)((context->adjacency[row * BLOCK_COUNT + word_index] >>
                      bit_index) &
                     UINT32_C(1));
}

static int partition_contains(const struct partition_vertex partition[],
                              int last, int vertex)
{
        int index;

        for (index = 0; index <= last; ++index) {
                if (vertex == partition[index].vertex) {
                        return 1;
                }
        }

        return 0;
}

static void apply_wqh(struct search_context *context,
                      const int neighbour_counts[])
{
        int vertex;
        int part_size = context->part_size;

        for (vertex = 0; vertex < context->vertex_count; ++vertex) {
                int row_offset;
                int word_index;
                int bit_index;
                int index;

                if (partition_contains(context->partition, 2 * part_size - 1,
                                       vertex)) {
                        continue;
                }

                if (!((neighbour_counts[2 * vertex] == part_size &&
                       neighbour_counts[2 * vertex + 1] == 0) ||
                      (neighbour_counts[2 * vertex] == 0 &&
                       neighbour_counts[2 * vertex + 1] == part_size))) {
                        continue;
                }

                row_offset = vertex * BLOCK_COUNT;
                word_index = vertex >> BLOCK_SHIFT;
                bit_index = vertex & BLOCK_MASK;

                for (index = 0; index < part_size; ++index) {
                        const struct partition_vertex *first =
                                &context->partition[index];
                        const struct partition_vertex *second =
                                &context->partition[index + part_size];
                        word_t vertex_bit = UINT32_C(1) << bit_index;

                        context->adjacency[row_offset + first->word_index] ^=
                                UINT32_C(1) << first->bit_index;
                        context->adjacency[row_offset + second->word_index] ^=
                                UINT32_C(1) << second->bit_index;
                        context->adjacency[first->row_offset + word_index] ^=
                                vertex_bit;
                        context->adjacency[second->row_offset + word_index] ^=
                                vertex_bit;
                }
        }
}

static enum graphswitching_result choose_partition(
        struct search_context *context, int current)
{
        int start = 0;
        int part_size = context->part_size;
        struct partition_vertex *selected = &context->partition[current];

        if (current == 2 * part_size) {
                return apply_partition_global(context);
        }

        if (current > 0) {
                start = context->partition[current - 1].vertex + 1;
                if (current == part_size) {
                        start = context->partition[0].vertex + 1;
                }
        }

        for (selected->vertex = start;
             selected->vertex < context->vertex_count;
             ++selected->vertex) {
                enum graphswitching_result result;

                if (partition_contains(context->partition, current - 1,
                                       selected->vertex)) {
                        continue;
                }

                selected->row_offset = selected->vertex * BLOCK_COUNT;
                selected->word_index = selected->vertex >> BLOCK_SHIFT;
                selected->bit_index = selected->vertex & BLOCK_MASK;

                if (current <= part_size / 2 || current == part_size) {
                        /* There is not enough information to prune yet. */
                } else if (current < part_size) {
                        if (!test_regular_11(context, current)) {
                                continue;
                        }
                } else {
                        if (!test_regular_21(context, current) ||
                            !test_regular_12(context, current) ||
                            !test_regular_22(context, current)) {
                                continue;
                        }
                }

                result = choose_partition(context, current + 1);
                if (result != GRAPHSWITCHING_SUCCESS) {
                        return result;
                }
        }

        return GRAPHSWITCHING_SUCCESS;
}

static int test_regular_11(const struct search_context *context, int current)
{
        int degrees[GRAPHSWITCHING_MAX_PART_SIZE] = {0};
        int minimum;
        int maximum;
        int index;

        get_degrees_11(context, current, degrees);
        minimum = degrees[0];
        maximum = degrees[0];

        for (index = 1; index <= current; ++index) {
                if (degrees[index] < minimum) {
                        minimum = degrees[index];
                } else if (degrees[index] > maximum) {
                        maximum = degrees[index];
                }

                if (maximum - minimum >
                    context->part_size - current - 1) {
                        return 0;
                }
        }

        return 1;
}

static int test_regular_22(const struct search_context *context, int current)
{
        int required[GRAPHSWITCHING_MAX_PART_SIZE] = {0};
        int degrees[GRAPHSWITCHING_MAX_PART_SIZE] = {0};
        int part_size = context->part_size;
        int index;

        get_degrees_11(context, part_size - 1, required);
        get_degrees_22(context, current, degrees);

        for (index = 0; index <= current - part_size; ++index) {
                if (required[0] < degrees[index] ||
                    required[0] >
                            degrees[index] + 2 * part_size - current - 1) {
                        return 0;
                }
        }

        return 1;
}

static int test_regular_12(const struct search_context *context, int current)
{
        int degrees[GRAPHSWITCHING_MAX_PART_SIZE] = {0};
        int minimum;
        int maximum;
        int part_size = context->part_size;
        int index;

        get_degrees_12(context, current, degrees);
        minimum = degrees[0];
        maximum = degrees[0];

        for (index = 1; index < part_size; ++index) {
                if (degrees[index] < minimum) {
                        minimum = degrees[index];
                } else if (degrees[index] > maximum) {
                        maximum = degrees[index];
                }

                if (maximum - minimum >
                    2 * part_size - current - 1) {
                        return 0;
                }
        }

        return 1;
}

static int test_regular_21(const struct search_context *context, int current)
{
        int degrees[GRAPHSWITCHING_MAX_PART_SIZE] = {0};
        int index;

        get_degrees_21(context, current, degrees);

        for (index = 1; index <= current - context->part_size; ++index) {
                if (degrees[index] != degrees[index - 1]) {
                        return 0;
                }
        }

        return 1;
}

static void get_degrees_11(const struct search_context *context, int current,
                           int degrees[])
{
        int part_size = context->part_size;
        int row_index;

        for (row_index = 0;
             row_index <= current && row_index < part_size;
             ++row_index) {
                int column_index;

                degrees[row_index] = 0;
                for (column_index = 0;
                     column_index <= current && column_index < part_size;
                     ++column_index) {
                        degrees[row_index] += adjacency_bit(
                                context,
                                context->partition[row_index].vertex,
                                context->partition[column_index].vertex);
                }
        }
}

static void get_degrees_22(const struct search_context *context, int current,
                           int degrees[])
{
        int part_size = context->part_size;
        int row_index;

        for (row_index = part_size; row_index <= current; ++row_index) {
                int column_index;

                degrees[row_index - part_size] = 0;
                for (column_index = part_size; column_index <= current;
                     ++column_index) {
                        degrees[row_index - part_size] += adjacency_bit(
                                context,
                                context->partition[row_index].vertex,
                                context->partition[column_index].vertex);
                }
        }
}

static void get_degrees_12(const struct search_context *context, int current,
                           int degrees[])
{
        int part_size = context->part_size;
        int row_index;

        for (row_index = 0; row_index < part_size; ++row_index) {
                int column_index;

                degrees[row_index] = 0;
                for (column_index = part_size; column_index <= current;
                     ++column_index) {
                        degrees[row_index] += adjacency_bit(
                                context,
                                context->partition[row_index].vertex,
                                context->partition[column_index].vertex);
                }
        }
}

static void get_degrees_21(const struct search_context *context, int current,
                           int degrees[])
{
        int part_size = context->part_size;
        int row_index;

        for (row_index = part_size; row_index <= current; ++row_index) {
                int column_index;

                degrees[row_index - part_size] = 0;
                for (column_index = 0; column_index < part_size;
                     ++column_index) {
                        degrees[row_index - part_size] += adjacency_bit(
                                context,
                                context->partition[row_index].vertex,
                                context->partition[column_index].vertex);
                }
        }
}

static enum graphswitching_result apply_partition_global(
        struct search_context *context)
{
        int neighbour_counts[2 * GRAPHSWITCHING_MAX_VERTICES];
        int any_switch = 0;
        int part_size = context->part_size;
        int vertex;

        for (vertex = 0; vertex < context->vertex_count; ++vertex) {
                int index;

                if (partition_contains(context->partition, 2 * part_size - 1,
                                       vertex)) {
                        continue;
                }

                neighbour_counts[2 * vertex] = 0;
                neighbour_counts[2 * vertex + 1] = 0;
                for (index = 0; index < part_size; ++index) {
                        neighbour_counts[2 * vertex] += adjacency_bit(
                                context, vertex,
                                context->partition[index].vertex);
                        neighbour_counts[2 * vertex + 1] += adjacency_bit(
                                context, vertex,
                                context->partition[index + part_size].vertex);
                }

                if (neighbour_counts[2 * vertex] ==
                    neighbour_counts[2 * vertex + 1]) {
                        continue;
                }

                if ((neighbour_counts[2 * vertex] == part_size &&
                     neighbour_counts[2 * vertex + 1] == 0) ||
                    (neighbour_counts[2 * vertex] == 0 &&
                     neighbour_counts[2 * vertex + 1] == part_size)) {
                        any_switch = 1;
                        continue;
                }

                return GRAPHSWITCHING_SUCCESS;
        }

        if (!any_switch) {
                return GRAPHSWITCHING_SUCCESS;
        }

        apply_wqh(context, neighbour_counts);
        {
                enum graphswitching_result result =
                        write_adjacency_matrix(context);

                apply_wqh(context, neighbour_counts);
                return result;
        }
}

static enum graphswitching_result write_adjacency_matrix(
        const struct search_context *context)
{
        char row[GRAPHSWITCHING_MAX_VERTICES + 1];
        int row_index;

        for (row_index = 0; row_index < context->vertex_count; ++row_index) {
                int column_index;

                for (column_index = 0;
                     column_index < context->vertex_count;
                     ++column_index) {
                        row[column_index] =
                                adjacency_bit(context, row_index, column_index)
                                        ? '1'
                                        : '0';
                }
                row[context->vertex_count] = '\n';

                if (fwrite(row, 1, (size_t)context->vertex_count + 1,
                           context->output) !=
                    (size_t)context->vertex_count + 1) {
                        return GRAPHSWITCHING_OUTPUT_ERROR;
                }
        }

        if (fputc('\n', context->output) == EOF) {
                return GRAPHSWITCHING_OUTPUT_ERROR;
        }

        return GRAPHSWITCHING_SUCCESS;
}
