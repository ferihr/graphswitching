/* Fixed-catalogue embedding search and switching. */

#include "graphswitching_internal.h"

#include <stdlib.h>
#include <string.h>

static enum graphswitching_result apply_fixed_method(
        struct graphswitching_search_context *context,
        const int selected[]);

static int subgraph_adjacency(uint64_t subgraph, int order,
                              int first, int second)
{
        int bit;

        if (first == second) {
                return 0;
        }
        if (first > second) {
                int temporary = first;

                first = second;
                second = temporary;
        }
        bit = first * (2 * order - first - 1) / 2 +
              second - first - 1;
        return (int)((subgraph >> bit) & UINT64_C(1));
}

static int build_method_blocks(
        const struct switching_method_definition *definition,
        uint16_t blocks[], uint16_t images[])
{
        int order = definition->order;
        int denominator = definition->denominator;
        int block_count = 0;
        unsigned int mask;

        for (mask = 0; mask < (1U << order); ++mask) {
                unsigned int image = 0;
                int valid = 1;
                int column;

                for (column = 0; column < order; ++column) {
                        int sum = 0;
                        int row;

                        for (row = 0; row < order; ++row) {
                                if ((mask >> row) & 1U) {
                                        sum += definition->numerator[
                                                row * order + column];
                                }
                        }
                        if (sum == denominator) {
                                image |= 1U << column;
                        } else if (sum != 0) {
                                valid = 0;
                                break;
                        }
                }
                if (valid) {
                        blocks[block_count] = (uint16_t)mask;
                        if (images != NULL) {
                                images[block_count] = (uint16_t)image;
                        }
                        ++block_count;
                }
        }
        return block_count;
}

static int fixed_method_pair_compatible(
        const struct graphswitching_search_context *context, uint64_t subgraph,
        int current, int graph_vertex)
{
        int previous;

        for (previous = 0; previous < current; ++previous) {
                if (graphswitching_adjacency_bit(
                            context, graph_vertex,
                            context->method_vertices[previous]) !=
                    subgraph_adjacency(
                            subgraph, context->definition->order,
                            current, previous)) {
                        return 0;
                }
        }
        return 1;
}

static int fixed_method_blocks_can_extend(
        const struct graphswitching_search_context *context, const int selected[],
        int selected_count, const uint16_t blocks[], int block_count)
{
        int order = context->definition->order;
        unsigned int known =
                selected_count == 0
                        ? 0
                        : (1U << selected_count) - 1U;
        int remaining = order - selected_count;
        int forced = 0;
        int vertex;

        for (vertex = 0; vertex < context->vertex_count; ++vertex) {
                unsigned int observed = 0;
                int is_selected = 0;
                int index;
                int block;

                for (index = 0; index < selected_count; ++index) {
                        if (selected[index] == vertex) {
                                is_selected = 1;
                                break;
                        }
                        if (graphswitching_adjacency_bit(
                                    context, vertex, selected[index])) {
                                observed |= 1U << index;
                        }
                }
                if (is_selected) {
                        continue;
                }
                for (block = 0; block < block_count; ++block) {
                        if (((unsigned int)blocks[block] & known) ==
                            observed) {
                                break;
                        }
                }
                if (block == block_count && ++forced > remaining) {
                        return 0;
                }
        }
        return 1;
}

static enum graphswitching_result choose_fixed_method_recursively(
        struct graphswitching_search_context *context,
        uint64_t subgraph, int current)
{
        int order = context->definition->order;
        int vertex;

        if (current == order) {
                return apply_fixed_method(
                        context, context->method_vertices);
        }

        for (vertex = 0; vertex < context->vertex_count; ++vertex) {
                enum graphswitching_result result;
                int index;

                for (index = 0; index < current; ++index) {
                        if (context->method_vertices[index] == vertex) {
                                break;
                        }
                }
                if (index < current ||
                    !fixed_method_pair_compatible(
                            context, subgraph, current, vertex)) {
                        continue;
                }

                context->method_vertices[current] = vertex;
                if (!fixed_method_blocks_can_extend(
                            context, context->method_vertices,
                            current + 1, context->method_blocks,
                            context->method_block_count)) {
                        continue;
                }
                result = choose_fixed_method_recursively(
                        context, subgraph, current + 1);
                if (result != GRAPHSWITCHING_SUCCESS) {
                        return result;
                }
        }
        return GRAPHSWITCHING_SUCCESS;
}

static enum graphswitching_result
choose_fixed_method_recursively_with_symmetry(
        struct graphswitching_search_context *context, uint64_t subgraph, int current)
{
        int order = context->definition->order;
        int symmetry_orbits[GRAPHSWITCHING_MAX_VERTICES];
        int vertex;

        if (current == order) {
                return apply_fixed_method(
                        context, context->method_vertices);
        }

        graphswitching_symmetry_fixed_orbits(
                context, current, symmetry_orbits);

        for (vertex = 0; vertex < context->vertex_count; ++vertex) {
                enum graphswitching_result result;
                int index;

                if (symmetry_orbits[vertex] != vertex) {
                        continue;
                }
                for (index = 0; index < current; ++index) {
                        if (context->method_vertices[index] == vertex) {
                                break;
                        }
                }
                if (index < current ||
                    !fixed_method_pair_compatible(
                            context, subgraph, current, vertex)) {
                        continue;
                }

                context->method_vertices[current] = vertex;
                if (!fixed_method_blocks_can_extend(
                            context, context->method_vertices,
                            current + 1, context->method_blocks,
                            context->method_block_count)) {
                        continue;
                }
                result = choose_fixed_method_recursively_with_symmetry(
                        context, subgraph, current + 1);
                if (result != GRAPHSWITCHING_SUCCESS) {
                        return result;
                }
        }
        return GRAPHSWITCHING_SUCCESS;
}

enum graphswitching_result graphswitching_choose_fixed_method(
        struct graphswitching_search_context *context)
{
        size_t subgraph_index;

        context->method_block_count = build_method_blocks(
                context->definition, context->method_blocks,
                context->method_images);
        for (subgraph_index = 0;
             subgraph_index <
                     context->definition->irreducible_subgraph_count;
             ++subgraph_index) {
                enum graphswitching_result result;
                uint64_t subgraph = context->definition
                                            ->irreducible_subgraphs[
                                                    subgraph_index];

                result = context->use_symmetry
                                 ? choose_fixed_method_recursively_with_symmetry(
                                           context, subgraph, 0)
                                 : choose_fixed_method_recursively(
                                           context, subgraph, 0);
                if (result != GRAPHSWITCHING_SUCCESS) {
                        return result;
                }
        }
        return GRAPHSWITCHING_SUCCESS;
}

static enum graphswitching_result apply_fixed_method(
        struct graphswitching_search_context *context, const int selected[])
{
        const struct switching_method_definition *definition =
                context->definition;
        size_t word_count =
                (size_t)context->vertex_count * GRAPHSWITCHING_BLOCK_COUNT;
        graphswitching_word_t *saved = (graphswitching_word_t *)malloc(word_count * sizeof(graphswitching_word_t));
        int order = definition->order;
        int denominator = definition->denominator;
        int denominator_squared = denominator * denominator;
        int changed = 0;
        int first;
        enum graphswitching_result result = GRAPHSWITCHING_SUCCESS;

        if (saved == NULL) {
                return GRAPHSWITCHING_MEMORY_ERROR;
        }
        memcpy(saved, context->adjacency, word_count * sizeof(graphswitching_word_t));

        for (first = 0; first < order; ++first) {
                int second;

                for (second = 0; second < order; ++second) {
                        int scaled = 0;
                        int row;
                        int value;

                        for (row = 0; row < order; ++row) {
                                int column;

                                for (column = 0; column < order; ++column) {
                                        int adjacent =
                                                (int)((saved[
                                                        selected[row] *
                                                                GRAPHSWITCHING_BLOCK_COUNT +
                                                        (selected[column] >>
                                                         GRAPHSWITCHING_BLOCK_SHIFT)] >>
                                                       (selected[column] &
                                                        GRAPHSWITCHING_BLOCK_MASK)) &
                                                      UINT32_C(1));

                                        scaled +=
                                                definition->numerator[
                                                        row * order + first] *
                                                adjacent *
                                                definition->numerator[
                                                        column * order +
                                                        second];
                                }
                        }
                        if (scaled != 0 &&
                            scaled != denominator_squared) {
                                result = GRAPHSWITCHING_INVALID_ARGUMENT;
                                goto restore;
                        }
                        value = scaled == denominator_squared;
                        if (value != graphswitching_adjacency_bit(
                                             context, selected[first],
                                             selected[second])) {
                                changed = 1;
                        }
                        graphswitching_set_adjacency_bit(
                                context, selected[first],
                                selected[second], value);
                }
        }

        for (first = 0; first < context->vertex_count; ++first) {
                int is_selected = 0;
                int index;

                for (index = 0; index < order; ++index) {
                        if (selected[index] == first) {
                                is_selected = 1;
                                break;
                        }
                }
                if (is_selected) {
                        continue;
                }
                for (index = 0; index < order; ++index) {
                        int scaled = 0;
                        int row;
                        int value;

                        for (row = 0; row < order; ++row) {
                                int adjacent =
                                        (int)((saved[
                                                first * GRAPHSWITCHING_BLOCK_COUNT +
                                                (selected[row] >>
                                                 GRAPHSWITCHING_BLOCK_SHIFT)] >>
                                               (selected[row] & GRAPHSWITCHING_BLOCK_MASK)) &
                                              UINT32_C(1));

                                scaled += definition->numerator[
                                                  row * order + index] *
                                          adjacent;
                        }
                        if (scaled != 0 && scaled != denominator) {
                                result = GRAPHSWITCHING_INVALID_ARGUMENT;
                                goto restore;
                        }
                        value = scaled == denominator;
                        if (value !=
                            (int)((saved[
                                    first * GRAPHSWITCHING_BLOCK_COUNT +
                                    (selected[index] >> GRAPHSWITCHING_BLOCK_SHIFT)] >>
                                   (selected[index] & GRAPHSWITCHING_BLOCK_MASK)) &
                                  UINT32_C(1))) {
                                changed = 1;
                        }
                        graphswitching_set_adjacency_bit(
                                context, first, selected[index], value);
                        graphswitching_set_adjacency_bit(
                                context, selected[index], first, value);
                }
        }

        if (changed) {
                result = graphswitching_write_graph(context);
        }

restore:
        memcpy(context->adjacency, saved, word_count * sizeof(graphswitching_word_t));
        free(saved);
        return result;
}
