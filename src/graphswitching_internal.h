/* Internal shared representation for graphswitching search modules. */

#ifndef GRAPHSWITCHING_INTERNAL_H
#define GRAPHSWITCHING_INTERNAL_H

#include "graphswitching.h"
#include "switching_methods.h"

#include <stdint.h>

#define GRAPHSWITCHING_WORD_BITS 32
#define GRAPHSWITCHING_BLOCK_SHIFT 5
#define GRAPHSWITCHING_BLOCK_MASK (GRAPHSWITCHING_WORD_BITS - 1)
/*
 * Keep the established row stride while enforcing the smaller public input
 * limit. Changing this layout is a separate performance change and makes the
 * fixed matcher's hottest address calculation slightly slower on this CPU.
 */
#define GRAPHSWITCHING_STORAGE_MAX_VERTICES 1536
#define GRAPHSWITCHING_BLOCK_COUNT \
        ((GRAPHSWITCHING_STORAGE_MAX_VERTICES + \
          GRAPHSWITCHING_WORD_BITS - 1) / \
         GRAPHSWITCHING_WORD_BITS)
#define GRAPHSWITCHING_MATRIX_BIT_CAPACITY \
        (GRAPHSWITCHING_MAX_VERTICES * GRAPHSWITCHING_MAX_VERTICES)
#define GRAPHSWITCHING_GM_SET_SIZE 4

typedef uint32_t graphswitching_word_t;

struct graphswitching_partition_vertex {
        int vertex;
        int row_offset;
        int word_index;
        int bit_index;
};

struct graphswitching_symmetry_search;

struct graphswitching_search_context {
        int vertex_count;
        int part_size;
        enum graphswitching_method method;
        const struct switching_method_definition *definition;
        graphswitching_word_t adjacency[
                GRAPHSWITCHING_MAX_VERTICES * GRAPHSWITCHING_BLOCK_COUNT];
        graphswitching_word_t c1_signatures[GRAPHSWITCHING_MAX_VERTICES];
        unsigned char c1_cross_degrees[GRAPHSWITCHING_MAX_VERTICES];
        int c1_degrees[GRAPHSWITCHING_MAX_PART_SIZE];
        int c2_degrees[GRAPHSWITCHING_MAX_PART_SIZE];
        int c1_to_c2_degrees[GRAPHSWITCHING_MAX_PART_SIZE];
        struct graphswitching_partition_vertex
                partition[2 * GRAPHSWITCHING_MAX_PART_SIZE];
        int method_vertices[GRAPHSWITCHING_MAX_METHOD_ORDER];
        uint16_t method_blocks[GRAPHSWITCHING_MAX_METHOD_BLOCKS];
        uint16_t method_images[GRAPHSWITCHING_MAX_METHOD_BLOCKS];
        int method_block_count;
        FILE *output;
        enum graphswitching_output_format output_format;
        enum graphswitching_symmetry_mode symmetry_mode;
        int use_symmetry;
        struct graphswitching_symmetry_search *symmetry;
};

static inline int graphswitching_adjacency_bit(
        const struct graphswitching_search_context *context,
        int row, int column)
{
        int word_index = column >> GRAPHSWITCHING_BLOCK_SHIFT;
        int bit_index = column & GRAPHSWITCHING_BLOCK_MASK;

        return (int)((context->adjacency[
                              row * GRAPHSWITCHING_BLOCK_COUNT + word_index] >>
                      bit_index) & UINT32_C(1));
}

static inline void graphswitching_set_adjacency_bit(
        struct graphswitching_search_context *context,
        int row, int column, int value)
{
        int word_index = column >> GRAPHSWITCHING_BLOCK_SHIFT;
        graphswitching_word_t bit =
                UINT32_C(1) << (column & GRAPHSWITCHING_BLOCK_MASK);
        graphswitching_word_t *entry = &context->adjacency[
                row * GRAPHSWITCHING_BLOCK_COUNT + word_index];

        if (value) {
                *entry |= bit;
        } else {
                *entry &= ~bit;
        }
}

enum graphswitching_result graphswitching_read_adjacency_matrix(
        struct graphswitching_search_context *context,
        FILE *input, int requested_vertices);
enum graphswitching_result graphswitching_write_output_header(
        const struct graphswitching_search_context *context);
enum graphswitching_result graphswitching_write_graph(
        const struct graphswitching_search_context *context);

enum graphswitching_result graphswitching_choose_gm_sets(
        struct graphswitching_search_context *context);
void graphswitching_build_gm_regular_completion_mask(
        const struct graphswitching_search_context *context,
        const int switching_set[GRAPHSWITCHING_GM_SET_SIZE],
        graphswitching_word_t mask[]);
enum graphswitching_result graphswitching_apply_gm_set(
        struct graphswitching_search_context *context,
        const int switching_set[GRAPHSWITCHING_GM_SET_SIZE]);

enum graphswitching_result graphswitching_choose_partition(
        struct graphswitching_search_context *context, int current);
enum graphswitching_result graphswitching_choose_fixed_method(
        struct graphswitching_search_context *context);

enum graphswitching_result graphswitching_symmetry_initialize(
        struct graphswitching_search_context *context);
void graphswitching_symmetry_destroy(
        struct graphswitching_search_context *context);
enum graphswitching_result graphswitching_choose_gm_sets_with_symmetry(
        struct graphswitching_search_context *context);
enum graphswitching_result graphswitching_symmetry_accept_gm_set(
        struct graphswitching_search_context *context,
        const int switching_set[GRAPHSWITCHING_GM_SET_SIZE],
        int *is_representative);
enum graphswitching_result graphswitching_symmetry_prepare_c1(
        struct graphswitching_search_context *context,
        int *is_representative);
enum graphswitching_result graphswitching_symmetry_prepare_c2(
        struct graphswitching_search_context *context,
        int *is_representative);
void graphswitching_symmetry_partition_orbits(
        struct graphswitching_search_context *context,
        int current, int orbits[]);
void graphswitching_symmetry_fixed_orbits(
        struct graphswitching_search_context *context,
        int current, int orbits[]);

#endif
