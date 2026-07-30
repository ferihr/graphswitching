/* Copyright (C) 2020--2023, Ferdinand Ihringer
 *
 * Generic Godsil--McKay and Wang--Qiu--Hu graph switching.
 *
 * The GM implementation follows code/gen_all_srgs.c. The WQH implementation
 * retains the algorithm of code/gen_all_srgs_wqh_generic.c, with cached C1
 * signatures, incremental degree state, and specialized bitset enumeration
 * for part sizes 3 through 5. Both use the same dynamic-order bit-packed
 * matrix and buffered row output.
 */

#include "graphswitching.h"

#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#define WORD_BITS 32
#define BLOCK_SHIFT 5
#define BLOCK_MASK (WORD_BITS - 1)
#define BLOCK_COUNT \
        ((GRAPHSWITCHING_MAX_VERTICES + WORD_BITS - 1) / WORD_BITS)
#define MATRIX_BIT_CAPACITY \
        (GRAPHSWITCHING_MAX_VERTICES * GRAPHSWITCHING_MAX_VERTICES)
#define GM_SET_SIZE 4

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
        word_t c1_signatures[GRAPHSWITCHING_MAX_VERTICES];
        unsigned char c1_cross_degrees[GRAPHSWITCHING_MAX_VERTICES];
        int c1_degrees[GRAPHSWITCHING_MAX_PART_SIZE];
        int c2_degrees[GRAPHSWITCHING_MAX_PART_SIZE];
        int c1_to_c2_degrees[GRAPHSWITCHING_MAX_PART_SIZE];
        struct partition_vertex
                partition[2 * GRAPHSWITCHING_MAX_PART_SIZE];
        FILE *output;
};

static int partition_contains(const struct partition_vertex partition[],
                              int last, int vertex);
static int gm_set_contains(const int switching_set[GM_SET_SIZE], int vertex);
static int adjacency_bit(const struct search_context *context, int row,
                         int column);
static void apply_gm(struct search_context *context,
                     const int switching_set[GM_SET_SIZE],
                     const int neighbour_counts[]);
static void apply_wqh(struct search_context *context,
                      const int neighbour_counts[]);
static enum graphswitching_result choose_gm_sets(
        struct search_context *context);
static enum graphswitching_result apply_gm_set(
        struct search_context *context,
        const int switching_set[GM_SET_SIZE]);
static enum graphswitching_result choose_partition(
        struct search_context *context, int current);
static enum graphswitching_result apply_partition_global(
        struct search_context *context);
static void cache_c1_signatures(struct search_context *context);
static void add_partition_vertex(struct search_context *context, int current);
static void remove_partition_vertex(struct search_context *context,
                                    int current);
static int build_specialized_candidate_mask(
        const struct search_context *context, int current, word_t mask[]);
static void intersect_with_neighbourhood(
        const struct search_context *context, word_t mask[], int vertex,
        int adjacent);
static void add_adjacency_pattern(
        const struct search_context *context, word_t destination[],
        int first_partition_index, int count, unsigned int pattern);
static int next_mask_vertex(const word_t mask[], int start, int limit);
static int test_regular_11(const struct search_context *context, int current);
static int test_regular_12(const struct search_context *context, int current);
static int test_regular_21(const struct search_context *context, int current);
static int test_regular_22(const struct search_context *context, int current);
static enum graphswitching_result read_adjacency_matrix(
        struct search_context *context, FILE *input, int requested_vertices);
static int parse_order_header(const char *line, int *vertex_count);
static int line_is_blank(const char *line);
static enum graphswitching_result collect_matrix_bits(
        FILE *input, unsigned char bits[], size_t *bit_count,
        int *declared_vertices);
static enum graphswitching_result write_adjacency_matrix(
        const struct search_context *context);

void graphswitching_options_init(struct graphswitching_options *options)
{
        if (options == NULL) {
                return;
        }

        options->method = GRAPHSWITCHING_METHOD_GM;
        options->vertex_count = 0;
        options->part_size = 2;
}

enum graphswitching_result graphswitching_generate_with_options(
        FILE *input,
        FILE *output,
        const struct graphswitching_options *options)
{
        struct search_context context;
        enum graphswitching_result result;

        if (input == NULL || output == NULL || options == NULL ||
            (options->method != GRAPHSWITCHING_METHOD_GM &&
             options->method != GRAPHSWITCHING_METHOD_WQH) ||
            options->vertex_count < 0 ||
            options->vertex_count > GRAPHSWITCHING_MAX_VERTICES ||
            (options->method == GRAPHSWITCHING_METHOD_WQH &&
             (options->part_size <= 0 ||
              options->part_size > GRAPHSWITCHING_MAX_PART_SIZE))) {
                return GRAPHSWITCHING_INVALID_ARGUMENT;
        }

        memset(&context, 0, sizeof(context));
        context.part_size = options->part_size;
        context.output = output;

        result = read_adjacency_matrix(&context, input,
                                       options->vertex_count);
        if (result != GRAPHSWITCHING_SUCCESS) {
                return result;
        }

        if (options->method == GRAPHSWITCHING_METHOD_WQH &&
            2 * context.part_size > context.vertex_count) {
                return GRAPHSWITCHING_INVALID_ARGUMENT;
        }

        if (fprintf(output, "n=%d\n", context.vertex_count) < 0) {
                return GRAPHSWITCHING_OUTPUT_ERROR;
        }

        if (options->method == GRAPHSWITCHING_METHOD_GM) {
                return choose_gm_sets(&context);
        }

        return choose_partition(&context, 0);
}

enum graphswitching_result graphswitching_generate(
        FILE *input,
        FILE *output,
        int vertex_count,
        int part_size)
{
        struct graphswitching_options options;

        graphswitching_options_init(&options);
        options.method = GRAPHSWITCHING_METHOD_WQH;
        options.vertex_count = vertex_count;
        options.part_size = part_size;
        return graphswitching_generate_with_options(input, output, &options);
}

const char *graphswitching_result_string(enum graphswitching_result result)
{
        switch (result) {
        case GRAPHSWITCHING_SUCCESS:
                return "success";
        case GRAPHSWITCHING_INVALID_ARGUMENT:
                return "invalid arguments";
        case GRAPHSWITCHING_INPUT_ERROR:
                return "could not read a square adjacency matrix";
        case GRAPHSWITCHING_OUTPUT_ERROR:
                return "could not write output";
        default:
                return "unknown error";
        }
}

static enum graphswitching_result read_adjacency_matrix(
        struct search_context *context, FILE *input, int requested_vertices)
{
        unsigned char bits[MATRIX_BIT_CAPACITY];
        size_t bit_count = 0;
        int declared_vertices = 0;
        int vertex_count = requested_vertices;
        size_t index;
        enum graphswitching_result result;

        result = collect_matrix_bits(input, bits, &bit_count,
                                     &declared_vertices);
        if (result != GRAPHSWITCHING_SUCCESS) {
                return result;
        }

        if (vertex_count == 0) {
                if (declared_vertices != 0) {
                        vertex_count = declared_vertices;
                } else {
                        for (vertex_count = 1;
                             vertex_count <= GRAPHSWITCHING_MAX_VERTICES;
                             ++vertex_count) {
                                if ((size_t)vertex_count *
                                            (size_t)vertex_count ==
                                    bit_count) {
                                        break;
                                }
                        }
                        if (vertex_count > GRAPHSWITCHING_MAX_VERTICES) {
                                return GRAPHSWITCHING_INPUT_ERROR;
                        }
                }
        } else if (declared_vertices != 0 &&
                   declared_vertices != vertex_count) {
                return GRAPHSWITCHING_INPUT_ERROR;
        }

        if (vertex_count <= 0 ||
            vertex_count > GRAPHSWITCHING_MAX_VERTICES ||
            bit_count !=
                    (size_t)vertex_count * (size_t)vertex_count) {
                return GRAPHSWITCHING_INPUT_ERROR;
        }

        context->vertex_count = vertex_count;
        for (index = 0; index < bit_count; ++index) {
                int row;
                int column;
                int word_index;
                int bit_index;

                if (bits[index] == 0) {
                        continue;
                }

                row = (int)(index / (size_t)vertex_count);
                column = (int)(index % (size_t)vertex_count);
                word_index = column >> BLOCK_SHIFT;
                bit_index = column & BLOCK_MASK;
                context->adjacency[row * BLOCK_COUNT + word_index] |=
                        UINT32_C(1) << bit_index;
        }

        return GRAPHSWITCHING_SUCCESS;
}

static enum graphswitching_result collect_matrix_bits(
        FILE *input, unsigned char bits[], size_t *bit_count,
        int *declared_vertices)
{
        char line[1024];
        int found_content = 0;
        int character;

        while (fgets(line, sizeof(line), input) != NULL) {
                const char *cursor;

                if (!found_content && line_is_blank(line)) {
                        continue;
                }

                if (!found_content) {
                        found_content = 1;
                        if (parse_order_header(line, declared_vertices)) {
                                break;
                        }
                }

                for (cursor = line; *cursor != '\0'; ++cursor) {
                        if (*cursor == '0' || *cursor == '1') {
                                if (*bit_count >= MATRIX_BIT_CAPACITY) {
                                        return GRAPHSWITCHING_INPUT_ERROR;
                                }
                                bits[(*bit_count)++] =
                                        (unsigned char)(*cursor - '0');
                        }
                }
                break;
        }

        if (ferror(input)) {
                return GRAPHSWITCHING_INPUT_ERROR;
        }

        while ((character = fgetc(input)) != EOF) {
                if (character == '0' || character == '1') {
                        if (*bit_count >= MATRIX_BIT_CAPACITY) {
                                return GRAPHSWITCHING_INPUT_ERROR;
                        }
                        bits[(*bit_count)++] =
                                (unsigned char)(character - '0');
                }
        }

        if (ferror(input) || *bit_count == 0) {
                return GRAPHSWITCHING_INPUT_ERROR;
        }

        return GRAPHSWITCHING_SUCCESS;
}

static int line_is_blank(const char *line)
{
        while (*line != '\0') {
                if (!isspace((unsigned char)*line)) {
                        return 0;
                }
                ++line;
        }

        return 1;
}

static int parse_order_header(const char *line, int *vertex_count)
{
        char *end;
        long parsed;

        while (isspace((unsigned char)*line)) {
                ++line;
        }
        if (*line++ != 'n') {
                return 0;
        }
        while (isspace((unsigned char)*line)) {
                ++line;
        }
        if (*line++ != '=') {
                return 0;
        }
        while (isspace((unsigned char)*line)) {
                ++line;
        }

        errno = 0;
        end = NULL;
        parsed = strtol(line, &end, 10);
        if (errno == ERANGE || end == line || parsed <= 0 ||
            parsed > INT_MAX) {
                return 0;
        }
        while (isspace((unsigned char)*end)) {
                ++end;
        }
        if (*end != '\0') {
                return 0;
        }

        *vertex_count = (int)parsed;
        return 1;
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

static int gm_set_contains(const int switching_set[GM_SET_SIZE], int vertex)
{
        int index;

        for (index = 0; index < GM_SET_SIZE; ++index) {
                if (switching_set[index] == vertex) {
                        return 1;
                }
        }

        return 0;
}

static void apply_gm(struct search_context *context,
                     const int switching_set[GM_SET_SIZE],
                     const int neighbour_counts[])
{
        int vertex;

        for (vertex = 0; vertex < context->vertex_count; ++vertex) {
                int row_offset;
                int word_index;
                int bit_index;
                int index;
                word_t vertex_bit;

                if (neighbour_counts[vertex] != GM_SET_SIZE / 2 ||
                    gm_set_contains(switching_set, vertex)) {
                        continue;
                }

                row_offset = vertex * BLOCK_COUNT;
                word_index = vertex >> BLOCK_SHIFT;
                bit_index = vertex & BLOCK_MASK;
                vertex_bit = UINT32_C(1) << bit_index;

                for (index = 0; index < GM_SET_SIZE; ++index) {
                        int switched_vertex = switching_set[index];
                        int switched_word = switched_vertex >> BLOCK_SHIFT;
                        int switched_bit = switched_vertex & BLOCK_MASK;

                        context->adjacency[row_offset + switched_word] ^=
                                UINT32_C(1) << switched_bit;
                        context->adjacency[switched_vertex * BLOCK_COUNT +
                                           word_index] ^= vertex_bit;
                }
        }
}

static enum graphswitching_result choose_gm_sets(
        struct search_context *context)
{
        int switching_set[GM_SET_SIZE];

        for (switching_set[0] = 0;
             switching_set[0] < context->vertex_count;
             ++switching_set[0]) {
                for (switching_set[1] = switching_set[0] + 1;
                     switching_set[1] < context->vertex_count;
                     ++switching_set[1]) {
                        for (switching_set[2] = switching_set[1] + 1;
                             switching_set[2] < context->vertex_count;
                             ++switching_set[2]) {
                                for (switching_set[3] =
                                             switching_set[2] + 1;
                                     switching_set[3] <
                                             context->vertex_count;
                                     ++switching_set[3]) {
                                        enum graphswitching_result result =
                                                apply_gm_set(
                                                        context,
                                                        switching_set);

                                        if (result !=
                                            GRAPHSWITCHING_SUCCESS) {
                                                return result;
                                        }
                                }
                        }
                }
        }

        return GRAPHSWITCHING_SUCCESS;
}

static enum graphswitching_result apply_gm_set(
        struct search_context *context,
        const int switching_set[GM_SET_SIZE])
{
        int neighbour_counts[GRAPHSWITCHING_MAX_VERTICES] = {0};
        int required_degree = 0;
        int any_switch = 0;
        int index;
        int vertex;

        for (index = 0; index < GM_SET_SIZE; ++index) {
                required_degree += adjacency_bit(
                        context, switching_set[0], switching_set[index]);
        }

        for (index = 1; index < GM_SET_SIZE; ++index) {
                int degree = 0;
                int column;

                for (column = 0; column < GM_SET_SIZE; ++column) {
                        degree += adjacency_bit(
                                context, switching_set[index],
                                switching_set[column]);
                }
                if (degree != required_degree) {
                        return GRAPHSWITCHING_SUCCESS;
                }
        }

        for (vertex = 0; vertex < context->vertex_count; ++vertex) {
                int column;

                if (gm_set_contains(switching_set, vertex)) {
                        continue;
                }

                for (column = 0; column < GM_SET_SIZE; ++column) {
                        neighbour_counts[vertex] += adjacency_bit(
                                context, vertex, switching_set[column]);
                }

                if (neighbour_counts[vertex] == GM_SET_SIZE / 2) {
                        any_switch = 1;
                } else if (neighbour_counts[vertex] != 0 &&
                           neighbour_counts[vertex] != GM_SET_SIZE) {
                        return GRAPHSWITCHING_SUCCESS;
                }
        }

        if (!any_switch) {
                return GRAPHSWITCHING_SUCCESS;
        }

        apply_gm(context, switching_set, neighbour_counts);
        {
                enum graphswitching_result result =
                        write_adjacency_matrix(context);

                apply_gm(context, switching_set, neighbour_counts);
                return result;
        }
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
        word_t candidate_mask[BLOCK_COUNT];
        int use_candidate_mask;
        struct partition_vertex *selected = &context->partition[current];

        if (current == 2 * part_size) {
                return apply_partition_global(context);
        }

        if (current == part_size) {
                cache_c1_signatures(context);
        }

        if (current > 0) {
                start = context->partition[current - 1].vertex + 1;
                if (current == part_size) {
                        start = context->partition[0].vertex + 1;
                }
        }

        use_candidate_mask = build_specialized_candidate_mask(
                context, current, candidate_mask);
        for (selected->vertex =
                     use_candidate_mask
                             ? next_mask_vertex(candidate_mask, start,
                                                context->vertex_count)
                             : start;
             selected->vertex < context->vertex_count;
             selected->vertex =
                     use_candidate_mask
                             ? next_mask_vertex(candidate_mask,
                                                selected->vertex + 1,
                                                context->vertex_count)
                             : selected->vertex + 1) {
                enum graphswitching_result result;
                int valid = 1;

                if (partition_contains(context->partition, current - 1,
                                       selected->vertex)) {
                        continue;
                }
                if (part_size > 0 && current > part_size &&
                    context->c1_cross_degrees[selected->vertex] !=
                            context->c1_cross_degrees[
                                    context->partition[part_size].vertex]) {
                        continue;
                }

                selected->row_offset = selected->vertex * BLOCK_COUNT;
                selected->word_index = selected->vertex >> BLOCK_SHIFT;
                selected->bit_index = selected->vertex & BLOCK_MASK;
                add_partition_vertex(context, current);

                if (current <= part_size / 2 || current == part_size) {
                        /* There is not enough information to prune yet. */
                } else if (current < part_size) {
                        if (!test_regular_11(context, current)) {
                                valid = 0;
                        }
                } else {
                        if (!test_regular_21(context, current) ||
                            !test_regular_12(context, current) ||
                            !test_regular_22(context, current)) {
                                valid = 0;
                        }
                }

                if (valid) {
                        result = choose_partition(context, current + 1);
                } else {
                        result = GRAPHSWITCHING_SUCCESS;
                }
                remove_partition_vertex(context, current);
                if (result != GRAPHSWITCHING_SUCCESS) {
                        return result;
                }
        }

        return GRAPHSWITCHING_SUCCESS;
}

static void cache_c1_signatures(struct search_context *context)
{
        int vertex;

        for (vertex = 0; vertex < context->vertex_count; ++vertex) {
                word_t signature = 0;
                int degree = 0;
                int index;

                for (index = 0; index < context->part_size; ++index) {
                        if (adjacency_bit(
                                    context, vertex,
                                    context->partition[index].vertex)) {
                                signature |= UINT32_C(1) << index;
                                ++degree;
                        }
                }

                context->c1_signatures[vertex] = signature;
                context->c1_cross_degrees[vertex] =
                        (unsigned char)degree;
        }
}

static void add_partition_vertex(struct search_context *context, int current)
{
        int part_size = context->part_size;
        int index;
        word_t signature;

        if (current < part_size) {
                context->c1_degrees[current] = 0;
                for (index = 0; index < current; ++index) {
                        int adjacent = adjacency_bit(
                                context,
                                context->partition[current].vertex,
                                context->partition[index].vertex);

                        context->c1_degrees[current] += adjacent;
                        context->c1_degrees[index] += adjacent;
                }
                return;
        }

        current -= part_size;
        context->c2_degrees[current] = 0;
        signature = context->c1_signatures[
                context->partition[current + part_size].vertex];
        for (index = 0; index < part_size; ++index) {
                context->c1_to_c2_degrees[index] +=
                        (int)((signature >> index) & UINT32_C(1));
        }
        for (index = 0; index < current; ++index) {
                int adjacent = adjacency_bit(
                        context,
                        context->partition[current + part_size].vertex,
                        context->partition[index + part_size].vertex);

                context->c2_degrees[current] += adjacent;
                context->c2_degrees[index] += adjacent;
        }
}

static void remove_partition_vertex(struct search_context *context,
                                    int current)
{
        int part_size = context->part_size;
        int index;
        word_t signature;

        if (current < part_size) {
                for (index = 0; index < current; ++index) {
                        int adjacent = adjacency_bit(
                                context,
                                context->partition[current].vertex,
                                context->partition[index].vertex);

                        context->c1_degrees[index] -= adjacent;
                }
                return;
        }

        current -= part_size;
        signature = context->c1_signatures[
                context->partition[current + part_size].vertex];
        for (index = 0; index < part_size; ++index) {
                context->c1_to_c2_degrees[index] -=
                        (int)((signature >> index) & UINT32_C(1));
        }
        for (index = 0; index < current; ++index) {
                int adjacent = adjacency_bit(
                        context,
                        context->partition[current + part_size].vertex,
                        context->partition[index + part_size].vertex);

                context->c2_degrees[index] -= adjacent;
        }
}

static int build_specialized_candidate_mask(
        const struct search_context *context, int current, word_t mask[])
{
        int part_size = context->part_size;
        int word_index;

        if (part_size < 3 || part_size > 5) {
                return 0;
        }

        for (word_index = 0; word_index < BLOCK_COUNT; ++word_index) {
                mask[word_index] = UINT32_MAX;
        }

        if (current == part_size - 1) {
                word_t allowed[BLOCK_COUNT] = {0};
                int selected_count = part_size - 1;
                unsigned int pattern;

                for (pattern = 0;
                     pattern < (UINT32_C(1) << selected_count);
                     ++pattern) {
                        int required_degree = 0;
                        int valid = 1;
                        int index;

                        for (index = 0; index < selected_count; ++index) {
                                required_degree +=
                                        (int)((pattern >> index) & 1U);
                        }
                        for (index = 0; index < selected_count; ++index) {
                                int adjacent =
                                        (int)((pattern >> index) & 1U);

                                if (context->c1_degrees[index] + adjacent !=
                                    required_degree) {
                                        valid = 0;
                                        break;
                                }
                        }
                        if (valid) {
                                add_adjacency_pattern(
                                        context, allowed, 0, selected_count,
                                        pattern);
                        }
                }

                memcpy(mask, allowed, sizeof(allowed));
                return 1;
        }

        if (current > part_size && current < 2 * part_size) {
                word_t allowed[BLOCK_COUNT] = {0};
                int selected_count = current - part_size;
                int required_degree = context->c1_degrees[0];
                unsigned int pattern;

                for (pattern = 0;
                     pattern < (UINT32_C(1) << selected_count);
                     ++pattern) {
                        int remaining = part_size - selected_count - 1;
                        int candidate_degree = 0;
                        int valid = 1;
                        int index;

                        for (index = 0; index < selected_count; ++index) {
                                int adjacent =
                                        (int)((pattern >> index) & 1U);
                                int degree =
                                        context->c2_degrees[index] +
                                        adjacent;

                                candidate_degree += adjacent;
                                if (degree > required_degree ||
                                    degree + remaining < required_degree) {
                                        valid = 0;
                                        break;
                                }
                        }
                        if (!valid ||
                            candidate_degree > required_degree ||
                            candidate_degree + remaining <
                                    required_degree) {
                                continue;
                        }

                        add_adjacency_pattern(
                                context, allowed, part_size,
                                selected_count, pattern);
                }

                memcpy(mask, allowed, sizeof(allowed));
                return 1;
        }

        return 0;
}

static void intersect_with_neighbourhood(
        const struct search_context *context, word_t mask[], int vertex,
        int adjacent)
{
        int row_offset = vertex * BLOCK_COUNT;
        int word_index;

        for (word_index = 0; word_index < BLOCK_COUNT; ++word_index) {
                word_t neighbours =
                        context->adjacency[row_offset + word_index];

                mask[word_index] &=
                        adjacent ? neighbours : ~neighbours;
        }
}

static void add_adjacency_pattern(
        const struct search_context *context, word_t destination[],
        int first_partition_index, int count, unsigned int pattern)
{
        word_t matching[BLOCK_COUNT];
        int index;
        int word_index;

        for (word_index = 0; word_index < BLOCK_COUNT; ++word_index) {
                matching[word_index] = UINT32_MAX;
        }
        for (index = 0; index < count; ++index) {
                intersect_with_neighbourhood(
                        context, matching,
                        context->partition[
                                first_partition_index + index].vertex,
                        (int)((pattern >> index) & 1U));
        }
        for (word_index = 0; word_index < BLOCK_COUNT; ++word_index) {
                destination[word_index] |= matching[word_index];
        }
}

static int next_mask_vertex(const word_t mask[], int start, int limit)
{
        int vertex;

        for (vertex = start; vertex < limit; ++vertex) {
                if (((mask[vertex >> BLOCK_SHIFT] >>
                      (vertex & BLOCK_MASK)) &
                     UINT32_C(1)) != 0) {
                        return vertex;
                }
        }

        return limit;
}

static int test_regular_11(const struct search_context *context, int current)
{
        int minimum;
        int maximum;
        int index;

        minimum = context->c1_degrees[0];
        maximum = context->c1_degrees[0];

        for (index = 1; index <= current; ++index) {
                if (context->c1_degrees[index] < minimum) {
                        minimum = context->c1_degrees[index];
                } else if (context->c1_degrees[index] > maximum) {
                        maximum = context->c1_degrees[index];
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
        int part_size = context->part_size;
        int index;

        for (index = 0; index <= current - part_size; ++index) {
                if (context->c1_degrees[0] <
                            context->c2_degrees[index] ||
                    context->c1_degrees[0] >
                            context->c2_degrees[index] +
                                    2 * part_size - current - 1) {
                        return 0;
                }
        }

        return 1;
}

static int test_regular_12(const struct search_context *context, int current)
{
        int minimum;
        int maximum;
        int part_size = context->part_size;
        int index;

        minimum = context->c1_to_c2_degrees[0];
        maximum = context->c1_to_c2_degrees[0];

        for (index = 1; index < part_size; ++index) {
                if (context->c1_to_c2_degrees[index] < minimum) {
                        minimum = context->c1_to_c2_degrees[index];
                } else if (context->c1_to_c2_degrees[index] > maximum) {
                        maximum = context->c1_to_c2_degrees[index];
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
        int part_size = context->part_size;
        int selected = current - part_size;

        return context->c1_cross_degrees[
                       context->partition[selected + part_size].vertex] ==
               context->c1_cross_degrees[
                       context->partition[part_size].vertex];
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

                neighbour_counts[2 * vertex] =
                        context->c1_cross_degrees[vertex];
                neighbour_counts[2 * vertex + 1] = 0;
                for (index = 0; index < part_size; ++index) {
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
