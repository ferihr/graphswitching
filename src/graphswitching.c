/* Copyright (C) 2020--2023, Ferdinand Ihringer
 *
 * Generic Godsil--McKay and Wang--Qiu--Hu graph switching.
 *
 * The GM implementation follows code/gen_all_srgs.c. The WQH implementation
 * retains the algorithm of code/gen_all_srgs_wqh_generic.c, with cached C1
 * signatures, incremental degree state, and specialized bitset enumeration
 * for part sizes 3 through 5. Both use the same dynamic-order bit-packed
 * matrix and buffered row output. A nauty-enabled build can additionally
 * restrict WQH enumeration to automorphism-orbit representatives.
 */

#include "graphswitching.h"

#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#ifdef GRAPHSWITCHING_WITH_NAUTY
#include <nauty.h>
#endif

#define WORD_BITS 32
#define BLOCK_SHIFT 5
#define BLOCK_MASK (WORD_BITS - 1)
#define BLOCK_COUNT \
        ((GRAPHSWITCHING_MAX_VERTICES + WORD_BITS - 1) / WORD_BITS)
#define MATRIX_BIT_CAPACITY \
        (GRAPHSWITCHING_MAX_VERTICES * GRAPHSWITCHING_MAX_VERTICES)
#define GM_SET_SIZE 4

typedef uint32_t word_t;

#ifdef GRAPHSWITCHING_WITH_NAUTY
struct permutation_group {
        int vertex_count;
        size_t generator_count;
        size_t generator_capacity;
        int *generators;
        int allocation_failed;
};

struct rank_set {
        uint64_t *slots;
        size_t count;
        size_t capacity;
};

struct subset_queue {
        int *vertices;
        size_t count;
        size_t capacity;
};

struct automorphism_search {
        struct permutation_group graph_group;
        struct permutation_group c1_stabilizer;
        struct rank_set seen_c1;
        struct rank_set seen_c2;
        struct subset_queue orbit_queue;
        uint64_t binomial[GRAPHSWITCHING_MAX_VERTICES + 1]
                         [GRAPHSWITCHING_MAX_PART_SIZE + 1];
};
#endif

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
        int use_symmetry;
#ifdef GRAPHSWITCHING_WITH_NAUTY
        struct automorphism_search automorphisms;
#endif
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
#ifdef GRAPHSWITCHING_WITH_NAUTY
static enum graphswitching_result initialize_automorphism_search(
        struct search_context *context);
static void destroy_automorphism_search(struct search_context *context);
static enum graphswitching_result prepare_c1_orbit(
        struct search_context *context, int *is_representative);
static enum graphswitching_result prepare_c2_orbit(
        struct search_context *context, int *is_representative);
static enum graphswitching_result compute_automorphism_group(
        const struct search_context *context,
        const struct partition_vertex fixed_subset[],
        int fixed_count,
        struct permutation_group *group);
static void destroy_permutation_group(struct permutation_group *group);
static enum graphswitching_result mark_subset_orbit(
        struct search_context *context,
        const struct permutation_group *group,
        int first_partition_index,
        struct rank_set *seen,
        int *is_representative);
static uint64_t subset_rank(const struct search_context *context,
                            const int subset[]);
static int rank_set_insert(struct rank_set *set, uint64_t rank);
static int rank_set_grow(struct rank_set *set);
static void rank_set_clear(struct rank_set *set);
static void rank_set_destroy(struct rank_set *set);
static int subset_queue_append(struct subset_queue *queue,
                               const int subset[], int part_size);
static void sort_subset(int subset[], int count);
#endif

void graphswitching_options_init(struct graphswitching_options *options)
{
        if (options == NULL) {
                return;
        }

        options->method = GRAPHSWITCHING_METHOD_GM;
        options->vertex_count = 0;
        options->part_size = 2;
        options->use_symmetry = 0;
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
            (options->use_symmetry != 0 &&
             options->use_symmetry != 1) ||
            (options->method != GRAPHSWITCHING_METHOD_WQH &&
             options->use_symmetry) ||
            (options->method == GRAPHSWITCHING_METHOD_WQH &&
             (options->part_size <= 0 ||
              options->part_size > GRAPHSWITCHING_MAX_PART_SIZE))) {
                return GRAPHSWITCHING_INVALID_ARGUMENT;
        }

#ifndef GRAPHSWITCHING_WITH_NAUTY
        if (options->use_symmetry) {
                return GRAPHSWITCHING_FEATURE_UNAVAILABLE;
        }
#endif

        memset(&context, 0, sizeof(context));
        context.part_size = options->part_size;
        context.output = output;
        context.use_symmetry = options->use_symmetry;

        result = read_adjacency_matrix(&context, input,
                                       options->vertex_count);
        if (result != GRAPHSWITCHING_SUCCESS) {
                return result;
        }

        if (options->method == GRAPHSWITCHING_METHOD_WQH &&
            2 * context.part_size > context.vertex_count) {
                return GRAPHSWITCHING_INVALID_ARGUMENT;
        }

#ifdef GRAPHSWITCHING_WITH_NAUTY
        if (context.use_symmetry) {
                result = initialize_automorphism_search(&context);
                if (result != GRAPHSWITCHING_SUCCESS) {
                        destroy_automorphism_search(&context);
                        return result;
                }
        }
#endif

        if (fprintf(output, "n=%d\n", context.vertex_count) < 0) {
#ifdef GRAPHSWITCHING_WITH_NAUTY
                destroy_automorphism_search(&context);
#endif
                return GRAPHSWITCHING_OUTPUT_ERROR;
        }

        if (options->method == GRAPHSWITCHING_METHOD_GM) {
                result = choose_gm_sets(&context);
        } else {
                result = choose_partition(&context, 0);
        }

#ifdef GRAPHSWITCHING_WITH_NAUTY
        destroy_automorphism_search(&context);
#endif
        return result;
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
        case GRAPHSWITCHING_FEATURE_UNAVAILABLE:
                return "symmetry search requires nauty support in this build";
        case GRAPHSWITCHING_MEMORY_ERROR:
                return "not enough memory";
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
#ifdef GRAPHSWITCHING_WITH_NAUTY
                if (context->use_symmetry) {
                        enum graphswitching_result result;
                        int is_representative;

                        result = prepare_c2_orbit(
                                context, &is_representative);
                        if (result != GRAPHSWITCHING_SUCCESS ||
                            !is_representative) {
                                return result;
                        }
                }
#endif
                return apply_partition_global(context);
        }

        if (current == part_size) {
#ifdef GRAPHSWITCHING_WITH_NAUTY
                if (context->use_symmetry) {
                        enum graphswitching_result result;
                        int is_representative;

                        result = prepare_c1_orbit(
                                context, &is_representative);
                        if (result != GRAPHSWITCHING_SUCCESS ||
                            !is_representative) {
                                return result;
                        }
                }
#endif
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

#ifdef GRAPHSWITCHING_WITH_NAUTY
static struct permutation_group *generator_target;

static void collect_automorphism_generator(
        int count, int *permutation, int *orbits, int numorbits,
        int stabilizer_vertex, int vertex_count)
{
        struct permutation_group *group = generator_target;
        int *resized;
        size_t new_capacity;

        (void)count;
        (void)orbits;
        (void)numorbits;
        (void)stabilizer_vertex;

        if (group == NULL || group->allocation_failed) {
                return;
        }

        if (group->generator_count == group->generator_capacity) {
                new_capacity = group->generator_capacity == 0
                                       ? 8
                                       : 2 * group->generator_capacity;
                resized = (int *)realloc(
                        group->generators,
                        new_capacity * (size_t)vertex_count * sizeof(int));
                if (resized == NULL) {
                        group->allocation_failed = 1;
                        return;
                }
                group->generators = resized;
                group->generator_capacity = new_capacity;
        }

        memcpy(group->generators +
                       group->generator_count * (size_t)vertex_count,
               permutation, (size_t)vertex_count * sizeof(int));
        ++group->generator_count;
}

static enum graphswitching_result compute_automorphism_group(
        const struct search_context *context,
        const struct partition_vertex fixed_subset[],
        int fixed_count,
        struct permutation_group *group)
{
        DEFAULTOPTIONS_GRAPH(options);
        statsblk statistics;
        graph *nauty_graph;
        int *lab;
        int *ptn;
        int *orbits;
        unsigned char fixed[GRAPHSWITCHING_MAX_VERTICES] = {0};
        int vertex_count = context->vertex_count;
        int set_words = SETWORDSNEEDED(vertex_count);
        int vertex;
        int position;

        destroy_permutation_group(group);
        group->vertex_count = vertex_count;

        nauty_graph = (graph *)calloc(
                (size_t)vertex_count * (size_t)set_words, sizeof(graph));
        lab = (int *)malloc((size_t)vertex_count * sizeof(int));
        ptn = (int *)malloc((size_t)vertex_count * sizeof(int));
        orbits = (int *)malloc((size_t)vertex_count * sizeof(int));
        if (nauty_graph == NULL || lab == NULL || ptn == NULL ||
            orbits == NULL) {
                free(nauty_graph);
                free(lab);
                free(ptn);
                free(orbits);
                return GRAPHSWITCHING_MEMORY_ERROR;
        }

        for (vertex = 0; vertex < vertex_count; ++vertex) {
                int neighbour;
                set *row = GRAPHROW(nauty_graph, vertex, set_words);

                for (neighbour = 0; neighbour < vertex_count; ++neighbour) {
                        if (adjacency_bit(context, vertex, neighbour)) {
                                ADDELEMENT(row, neighbour);
                        }
                }
        }

        if (fixed_count > 0) {
                position = 0;
                for (vertex = 0; vertex < fixed_count; ++vertex) {
                        int selected = fixed_subset[vertex].vertex;

                        fixed[selected] = 1;
                        lab[position] = selected;
                        ptn[position] =
                                vertex + 1 < fixed_count ? 1 : 0;
                        ++position;
                }
                for (vertex = 0; vertex < vertex_count; ++vertex) {
                        if (!fixed[vertex]) {
                                lab[position] = vertex;
                                ptn[position] =
                                        position + 1 < vertex_count ? 1 : 0;
                                ++position;
                        }
                }
                options.defaultptn = FALSE;
        }

        options.userautomproc = collect_automorphism_generator;
        options.schreier = TRUE;
        generator_target = group;
        nauty_check(WORDSIZE, set_words, vertex_count, NAUTYVERSIONID);
        densenauty(nauty_graph, lab, ptn, orbits, &options, &statistics,
                   set_words, vertex_count, NULL);
        generator_target = NULL;

        free(nauty_graph);
        free(lab);
        free(ptn);
        free(orbits);

        if (group->allocation_failed) {
                destroy_permutation_group(group);
                return GRAPHSWITCHING_MEMORY_ERROR;
        }
        if (statistics.errstatus != 0) {
                destroy_permutation_group(group);
                return GRAPHSWITCHING_INVALID_ARGUMENT;
        }

        return GRAPHSWITCHING_SUCCESS;
}

static void destroy_permutation_group(struct permutation_group *group)
{
        free(group->generators);
        memset(group, 0, sizeof(*group));
}

static enum graphswitching_result initialize_automorphism_search(
        struct search_context *context)
{
        struct automorphism_search *search = &context->automorphisms;
        int vertex;
        int part;

        for (vertex = 0;
             vertex <= GRAPHSWITCHING_MAX_VERTICES;
             ++vertex) {
                search->binomial[vertex][0] = UINT64_C(1);
                for (part = 1;
                     part <= GRAPHSWITCHING_MAX_PART_SIZE;
                     ++part) {
                        if (part > vertex) {
                                search->binomial[vertex][part] = 0;
                        } else if (part == vertex) {
                                search->binomial[vertex][part] =
                                        UINT64_C(1);
                        } else {
                                search->binomial[vertex][part] =
                                        search->binomial[vertex - 1][part] +
                                        search->binomial[vertex - 1]
                                                                [part - 1];
                        }
                }
        }

        return compute_automorphism_group(
                context, NULL, 0, &search->graph_group);
}

static void destroy_automorphism_search(struct search_context *context)
{
        struct automorphism_search *search = &context->automorphisms;

        destroy_permutation_group(&search->graph_group);
        destroy_permutation_group(&search->c1_stabilizer);
        rank_set_destroy(&search->seen_c1);
        rank_set_destroy(&search->seen_c2);
        free(search->orbit_queue.vertices);
        memset(&search->orbit_queue, 0, sizeof(search->orbit_queue));
}

static enum graphswitching_result prepare_c1_orbit(
        struct search_context *context, int *is_representative)
{
        struct automorphism_search *search = &context->automorphisms;
        enum graphswitching_result result;

        result = mark_subset_orbit(
                context, &search->graph_group, 0,
                &search->seen_c1, is_representative);
        if (result != GRAPHSWITCHING_SUCCESS || !*is_representative) {
                return result;
        }

        rank_set_clear(&search->seen_c2);
        return compute_automorphism_group(
                context, context->partition, context->part_size,
                &search->c1_stabilizer);
}

static enum graphswitching_result prepare_c2_orbit(
        struct search_context *context, int *is_representative)
{
        struct automorphism_search *search = &context->automorphisms;

        return mark_subset_orbit(
                context, &search->c1_stabilizer, context->part_size,
                &search->seen_c2, is_representative);
}

static enum graphswitching_result mark_subset_orbit(
        struct search_context *context,
        const struct permutation_group *group,
        int first_partition_index,
        struct rank_set *seen,
        int *is_representative)
{
        struct subset_queue *queue =
                &context->automorphisms.orbit_queue;
        int subset[GRAPHSWITCHING_MAX_PART_SIZE];
        int image[GRAPHSWITCHING_MAX_PART_SIZE];
        int part_size = context->part_size;
        int inserted;
        int index;
        size_t head;

        for (index = 0; index < part_size; ++index) {
                subset[index] =
                        context->partition[
                                first_partition_index + index].vertex;
        }

        inserted = rank_set_insert(seen, subset_rank(context, subset));
        if (inserted < 0) {
                return GRAPHSWITCHING_MEMORY_ERROR;
        }
        if (!inserted) {
                *is_representative = 0;
                return GRAPHSWITCHING_SUCCESS;
        }

        queue->count = 0;
        if (!subset_queue_append(queue, subset, part_size)) {
                return GRAPHSWITCHING_MEMORY_ERROR;
        }

        for (head = 0; head < queue->count; ++head) {
                size_t generator;

                memcpy(subset,
                       queue->vertices + head * (size_t)part_size,
                       (size_t)part_size * sizeof(int));
                for (generator = 0;
                     generator < group->generator_count;
                     ++generator) {
                        const int *permutation =
                                group->generators +
                                generator *
                                        (size_t)group->vertex_count;

                        for (index = 0; index < part_size; ++index) {
                                image[index] =
                                        permutation[subset[index]];
                        }
                        sort_subset(image, part_size);
                        inserted = rank_set_insert(
                                seen, subset_rank(context, image));
                        if (inserted < 0) {
                                return GRAPHSWITCHING_MEMORY_ERROR;
                        }
                        if (inserted &&
                            !subset_queue_append(
                                    queue, image, part_size)) {
                                return GRAPHSWITCHING_MEMORY_ERROR;
                        }
                }
        }

        *is_representative = 1;
        return GRAPHSWITCHING_SUCCESS;
}

static uint64_t subset_rank(const struct search_context *context,
                            const int subset[])
{
        uint64_t rank = 0;
        int index;

        for (index = 0; index < context->part_size; ++index) {
                rank += context->automorphisms
                                .binomial[subset[index]][index + 1];
        }
        return rank;
}

static uint64_t mix_rank(uint64_t rank)
{
        rank ^= rank >> 30;
        rank *= UINT64_C(0xbf58476d1ce4e5b9);
        rank ^= rank >> 27;
        rank *= UINT64_C(0x94d049bb133111eb);
        rank ^= rank >> 31;
        return rank;
}

static int rank_set_insert(struct rank_set *set, uint64_t rank)
{
        size_t slot;

        if (set->capacity == 0 ||
            10 * (set->count + 1) > 7 * set->capacity) {
                if (!rank_set_grow(set)) {
                        return -1;
                }
        }

        slot = (size_t)mix_rank(rank) & (set->capacity - 1);
        while (set->slots[slot] != UINT64_MAX) {
                if (set->slots[slot] == rank) {
                        return 0;
                }
                slot = (slot + 1) & (set->capacity - 1);
        }

        set->slots[slot] = rank;
        ++set->count;
        return 1;
}

static int rank_set_grow(struct rank_set *set)
{
        uint64_t *old_slots = set->slots;
        size_t old_capacity = set->capacity;
        size_t new_capacity =
                old_capacity == 0 ? 1024 : 2 * old_capacity;
        size_t index;

        set->slots = (uint64_t *)malloc(
                new_capacity * sizeof(uint64_t));
        if (set->slots == NULL) {
                set->slots = old_slots;
                return 0;
        }
        set->capacity = new_capacity;
        set->count = 0;
        memset(set->slots, 0xff,
               new_capacity * sizeof(uint64_t));

        for (index = 0; index < old_capacity; ++index) {
                if (old_slots[index] != UINT64_MAX) {
                        size_t slot =
                                (size_t)mix_rank(old_slots[index]) &
                                (new_capacity - 1);

                        while (set->slots[slot] != UINT64_MAX) {
                                slot = (slot + 1) &
                                       (new_capacity - 1);
                        }
                        set->slots[slot] = old_slots[index];
                        ++set->count;
                }
        }
        free(old_slots);
        return 1;
}

static void rank_set_clear(struct rank_set *set)
{
        if (set->slots != NULL) {
                memset(set->slots, 0xff,
                       set->capacity * sizeof(uint64_t));
        }
        set->count = 0;
}

static void rank_set_destroy(struct rank_set *set)
{
        free(set->slots);
        memset(set, 0, sizeof(*set));
}

static int subset_queue_append(struct subset_queue *queue,
                               const int subset[], int part_size)
{
        if (queue->count == queue->capacity) {
                size_t new_capacity =
                        queue->capacity == 0
                                ? 1024
                                : 2 * queue->capacity;
                int *resized;

                if (new_capacity >
                    SIZE_MAX / (size_t)part_size / sizeof(int)) {
                        return 0;
                }
                resized = (int *)realloc(
                        queue->vertices,
                        new_capacity * (size_t)part_size * sizeof(int));
                if (resized == NULL) {
                        return 0;
                }
                queue->vertices = resized;
                queue->capacity = new_capacity;
        }

        memcpy(queue->vertices +
                       queue->count * (size_t)part_size,
               subset, (size_t)part_size * sizeof(int));
        ++queue->count;
        return 1;
}

static void sort_subset(int subset[], int count)
{
        int index;

        for (index = 1; index < count; ++index) {
                int value = subset[index];
                int position = index;

                while (position > 0 &&
                       subset[position - 1] > value) {
                        subset[position] = subset[position - 1];
                        --position;
                }
                subset[position] = value;
        }
}
#endif

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
