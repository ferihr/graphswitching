/* Godsil--McKay four-vertex switching search. */

#include "graphswitching_internal.h"

#include <string.h>

static int gm_set_contains(const int switching_set[GRAPHSWITCHING_GM_SET_SIZE], int vertex)
{
        int index;

        for (index = 0; index < GRAPHSWITCHING_GM_SET_SIZE; ++index) {
                if (switching_set[index] == vertex) {
                        return 1;
                }
        }

        return 0;
}

static void apply_gm(struct graphswitching_search_context *context,
                     const int switching_set[GRAPHSWITCHING_GM_SET_SIZE],
                     const int neighbour_counts[])
{
        int vertex;

        for (vertex = 0; vertex < context->vertex_count; ++vertex) {
                int row_offset;
                int word_index;
                int bit_index;
                int index;
                graphswitching_word_t vertex_bit;

                if (neighbour_counts[vertex] != GRAPHSWITCHING_GM_SET_SIZE / 2 ||
                    gm_set_contains(switching_set, vertex)) {
                        continue;
                }

                row_offset = vertex * GRAPHSWITCHING_BLOCK_COUNT;
                word_index = vertex >> GRAPHSWITCHING_BLOCK_SHIFT;
                bit_index = vertex & GRAPHSWITCHING_BLOCK_MASK;
                vertex_bit = UINT32_C(1) << bit_index;

                for (index = 0; index < GRAPHSWITCHING_GM_SET_SIZE; ++index) {
                        int switched_vertex = switching_set[index];
                        int switched_word = switched_vertex >> GRAPHSWITCHING_BLOCK_SHIFT;
                        int switched_bit = switched_vertex & GRAPHSWITCHING_BLOCK_MASK;

                        context->adjacency[row_offset + switched_word] ^=
                                UINT32_C(1) << switched_bit;
                        context->adjacency[switched_vertex * GRAPHSWITCHING_BLOCK_COUNT +
                                           word_index] ^= vertex_bit;
                }
        }
}

enum graphswitching_result graphswitching_choose_gm_sets(
        struct graphswitching_search_context *context)
{
        int switching_set[GRAPHSWITCHING_GM_SET_SIZE];

        if (context->use_symmetry) {
                return graphswitching_choose_gm_sets_with_symmetry(context);
        }

        for (switching_set[0] = 0;
             switching_set[0] < context->vertex_count;
             ++switching_set[0]) {
                for (switching_set[1] = switching_set[0] + 1;
                     switching_set[1] < context->vertex_count;
                     ++switching_set[1]) {
                        for (switching_set[2] = switching_set[1] + 1;
                             switching_set[2] < context->vertex_count;
                             ++switching_set[2]) {
                                graphswitching_word_t completion_mask[GRAPHSWITCHING_BLOCK_COUNT];

                                graphswitching_build_gm_regular_completion_mask(
                                        context, switching_set,
                                        completion_mask);
                                for (switching_set[3] =
                                             switching_set[2] + 1;
                                     switching_set[3] <
                                             context->vertex_count;
                                     ++switching_set[3]) {
                                        if (((completion_mask[
                                                      switching_set[3] >>
                                                      GRAPHSWITCHING_BLOCK_SHIFT] >>
                                              (switching_set[3] &
                                               GRAPHSWITCHING_BLOCK_MASK)) &
                                             UINT32_C(1)) == 0) {
                                                continue;
                                        }
                                        enum graphswitching_result result =
                                                graphswitching_apply_gm_set(
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

void graphswitching_build_gm_regular_completion_mask(
        const struct graphswitching_search_context *context,
        const int switching_set[GRAPHSWITCHING_GM_SET_SIZE], graphswitching_word_t mask[])
{
        int degrees[GRAPHSWITCHING_GM_SET_SIZE - 1] = {0};
        int degree_sum = 0;
        int pattern = 0;
        int index;
        int other;
        int word_index;

        for (index = 0; index < GRAPHSWITCHING_GM_SET_SIZE - 1; ++index) {
                for (other = index + 1;
                     other < GRAPHSWITCHING_GM_SET_SIZE - 1;
                     ++other) {
                        int adjacent = graphswitching_adjacency_bit(
                                context, switching_set[index],
                                switching_set[other]);

                        degrees[index] += adjacent;
                        degrees[other] += adjacent;
                }
        }
        for (index = 0; index < GRAPHSWITCHING_GM_SET_SIZE - 1; ++index) {
                degree_sum += degrees[index];
        }

        if (degree_sum == 0) {
                pattern = 0;
        } else if (degree_sum == 6) {
                pattern = (1 << (GRAPHSWITCHING_GM_SET_SIZE - 1)) - 1;
        } else {
                for (index = 0; index < GRAPHSWITCHING_GM_SET_SIZE - 1; ++index) {
                        if ((degree_sum == 2 && degrees[index] == 0) ||
                            (degree_sum == 4 && degrees[index] == 1)) {
                                pattern |= 1 << index;
                        }
                }
        }

        for (word_index = 0; word_index < GRAPHSWITCHING_BLOCK_COUNT; ++word_index) {
                mask[word_index] = UINT32_MAX;
        }
        for (index = 0; index < GRAPHSWITCHING_GM_SET_SIZE - 1; ++index) {
                int row_offset = switching_set[index] * GRAPHSWITCHING_BLOCK_COUNT;

                for (word_index = 0;
                     word_index < GRAPHSWITCHING_BLOCK_COUNT;
                     ++word_index) {
                        graphswitching_word_t neighbours =
                                context->adjacency[
                                        row_offset + word_index];

                        mask[word_index] &=
                                ((pattern >> index) & 1)
                                        ? neighbours
                                        : ~neighbours;
                }
        }
}

enum graphswitching_result graphswitching_apply_gm_set(
        struct graphswitching_search_context *context,
        const int switching_set[GRAPHSWITCHING_GM_SET_SIZE])
{
        int neighbour_counts[GRAPHSWITCHING_MAX_VERTICES] = {0};
        int required_degree = 0;
        int any_switch = 0;
        int index;
        int vertex;

        for (index = 0; index < GRAPHSWITCHING_GM_SET_SIZE; ++index) {
                required_degree += graphswitching_adjacency_bit(
                        context, switching_set[0], switching_set[index]);
        }

        for (index = 1; index < GRAPHSWITCHING_GM_SET_SIZE; ++index) {
                int degree = 0;
                int column;

                for (column = 0; column < GRAPHSWITCHING_GM_SET_SIZE; ++column) {
                        degree += graphswitching_adjacency_bit(
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

                for (column = 0; column < GRAPHSWITCHING_GM_SET_SIZE; ++column) {
                        neighbour_counts[vertex] += graphswitching_adjacency_bit(
                                context, vertex, switching_set[column]);
                }

                if (neighbour_counts[vertex] == GRAPHSWITCHING_GM_SET_SIZE / 2) {
                        any_switch = 1;
                } else if (neighbour_counts[vertex] != 0 &&
                           neighbour_counts[vertex] != GRAPHSWITCHING_GM_SET_SIZE) {
                        return GRAPHSWITCHING_SUCCESS;
                }
        }

        if (!any_switch) {
                return GRAPHSWITCHING_SUCCESS;
        }

        if (context->use_symmetry) {
                enum graphswitching_result result;
                int is_representative;

                result = graphswitching_symmetry_accept_gm_set(
                        context, switching_set, &is_representative);
                if (result != GRAPHSWITCHING_SUCCESS ||
                    !is_representative) {
                        return result;
                }
        }

        apply_gm(context, switching_set, neighbour_counts);
        {
                enum graphswitching_result result =
                        graphswitching_write_graph(context);

                apply_gm(context, switching_set, neighbour_counts);
                return result;
        }
}
