/* Wang--Qiu--Hu partition search. */

#include "graphswitching_internal.h"

#include <string.h>

static void cache_c1_signatures(
        struct graphswitching_search_context *context);
static void add_partition_vertex(
        struct graphswitching_search_context *context, int current);
static void remove_partition_vertex(
        struct graphswitching_search_context *context, int current);
static int build_specialized_candidate_mask(
        const struct graphswitching_search_context *context,
        int current, graphswitching_word_t mask[]);
static void intersect_with_neighbourhood(
        const struct graphswitching_search_context *context,
        graphswitching_word_t mask[], int vertex, int adjacent);
static void add_adjacency_pattern(
        const struct graphswitching_search_context *context,
        graphswitching_word_t destination[],
        int first_partition_index, int count, unsigned int pattern);
static int next_mask_vertex(
        const graphswitching_word_t mask[], int start, int limit);
static int test_regular_11(
        const struct graphswitching_search_context *context, int current);
static int test_regular_12(
        const struct graphswitching_search_context *context, int current);
static int test_regular_21(
        const struct graphswitching_search_context *context, int current);
static int test_regular_22(
        const struct graphswitching_search_context *context, int current);
static int partition_blocks_can_extend(
        const struct graphswitching_search_context *context,
        int selected_count);
static enum graphswitching_result apply_partition_global(
        struct graphswitching_search_context *context);

static void apply_wqh(struct graphswitching_search_context *context,
                      const int neighbour_counts[])
{
        int vertex;
        int part_size = context->part_size;

        for (vertex = 0; vertex < context->vertex_count; ++vertex) {
                int row_offset;
                int word_index;
                int bit_index;
                int index;

                if (context->in_partition[vertex]) {
                        continue;
                }

                if (!((neighbour_counts[2 * vertex] == part_size &&
                       neighbour_counts[2 * vertex + 1] == 0) ||
                      (neighbour_counts[2 * vertex] == 0 &&
                       neighbour_counts[2 * vertex + 1] == part_size))) {
                        continue;
                }

                row_offset = vertex * GRAPHSWITCHING_BLOCK_COUNT;
                word_index = vertex >> GRAPHSWITCHING_BLOCK_SHIFT;
                bit_index = vertex & GRAPHSWITCHING_BLOCK_MASK;

                for (index = 0; index < part_size; ++index) {
                        const struct graphswitching_partition_vertex *first =
                                &context->partition[index];
                        const struct graphswitching_partition_vertex *second =
                                &context->partition[index + part_size];
                        graphswitching_word_t vertex_bit = UINT32_C(1) << bit_index;

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

enum graphswitching_result graphswitching_choose_partition(
        struct graphswitching_search_context *context, int current)
{
        int start = 0;
        int part_size = context->part_size;
        graphswitching_word_t candidate_mask[GRAPHSWITCHING_BLOCK_COUNT];
        int use_candidate_mask;
        struct graphswitching_partition_vertex *selected = &context->partition[current];
        int symmetry_orbits[GRAPHSWITCHING_MAX_VERTICES];

        if (current == 2 * part_size) {
                if (context->use_symmetry) {
                        enum graphswitching_result result;
                        int is_representative;

                        result = graphswitching_symmetry_prepare_c2(
                                context, &is_representative);
                        if (result != GRAPHSWITCHING_SUCCESS ||
                            !is_representative) {
                                return result;
                        }
                }
                return apply_partition_global(context);
        }

        if (current == part_size) {
                if (context->use_symmetry) {
                        enum graphswitching_result result;
                        int is_representative;

                        result = graphswitching_symmetry_prepare_c1(
                                context, &is_representative);
                        if (result != GRAPHSWITCHING_SUCCESS ||
                            !is_representative) {
                                return result;
                        }
                }
                cache_c1_signatures(context);
        }

        if (current > 0) {
                start = context->partition[current - 1].vertex + 1;
                if (current == part_size) {
                        start = context->partition[0].vertex + 1;
                }
        }

        if (context->use_symmetry) {
                start = 0;
                graphswitching_symmetry_partition_orbits(
                        context, current, symmetry_orbits);
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

                if (context->use_symmetry &&
                    symmetry_orbits[selected->vertex] !=
                            selected->vertex) {
                        continue;
                }
                if (context->in_partition[selected->vertex]) {
                        continue;
                }
                if (part_size > 0 && current > part_size &&
                    context->c1_cross_degrees[selected->vertex] !=
                            context->c1_cross_degrees[
                                    context->partition[part_size].vertex]) {
                        continue;
                }

                selected->row_offset = selected->vertex * GRAPHSWITCHING_BLOCK_COUNT;
                selected->word_index = selected->vertex >> GRAPHSWITCHING_BLOCK_SHIFT;
                selected->bit_index = selected->vertex & GRAPHSWITCHING_BLOCK_MASK;
                context->in_partition[selected->vertex] = 1;
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

                /*
                 * With at most one C2 vertex every observed outside degree
                 * can still extend to a balanced or extreme block. Delay
                 * this scan until it can reject a branch.
                 */
                if (valid && part_size >= 3 &&
                    (context->use_symmetry || part_size <= 4) &&
                    current + 1 > part_size + 1) {
                        if (!partition_blocks_can_extend(
                                    context, current + 1)) {
                                valid = 0;
                        }
                }

                if (valid) {
                        result = graphswitching_choose_partition(context, current + 1);
                } else {
                        result = GRAPHSWITCHING_SUCCESS;
                }
                remove_partition_vertex(context, current);
                context->in_partition[selected->vertex] = 0;
                if (result != GRAPHSWITCHING_SUCCESS) {
                        return result;
                }
        }

        return GRAPHSWITCHING_SUCCESS;
}
static void cache_c1_signatures(struct graphswitching_search_context *context)
{
        int vertex;

        for (vertex = 0; vertex < context->vertex_count; ++vertex) {
                graphswitching_word_t signature = 0;
                int degree = 0;
                int index;

                for (index = 0; index < context->part_size; ++index) {
                        if (graphswitching_adjacency_bit(
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

static void add_partition_vertex(struct graphswitching_search_context *context, int current)
{
        int part_size = context->part_size;
        int index;
        graphswitching_word_t signature;

        if (current < part_size) {
                context->c1_degrees[current] = 0;
                for (index = 0; index < current; ++index) {
                        int adjacent = graphswitching_adjacency_bit(
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
                int adjacent = graphswitching_adjacency_bit(
                        context,
                        context->partition[current + part_size].vertex,
                        context->partition[index + part_size].vertex);

                context->c2_degrees[current] += adjacent;
                context->c2_degrees[index] += adjacent;
        }
}

static void remove_partition_vertex(struct graphswitching_search_context *context,
                                    int current)
{
        int part_size = context->part_size;
        int index;
        graphswitching_word_t signature;

        if (current < part_size) {
                for (index = 0; index < current; ++index) {
                        int adjacent = graphswitching_adjacency_bit(
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
                int adjacent = graphswitching_adjacency_bit(
                        context,
                        context->partition[current + part_size].vertex,
                        context->partition[index + part_size].vertex);

                context->c2_degrees[index] -= adjacent;
        }
}

static int build_specialized_candidate_mask(
        const struct graphswitching_search_context *context, int current, graphswitching_word_t mask[])
{
        int part_size = context->part_size;
        int word_index;

        /* At p=6 a completion still has at most 2^(p-1) = 32 patterns. */
        if (part_size < 3 || part_size > 6) {
                return 0;
        }

        for (word_index = 0; word_index < GRAPHSWITCHING_BLOCK_COUNT; ++word_index) {
                mask[word_index] = UINT32_MAX;
        }

        if (current == part_size - 1) {
                graphswitching_word_t allowed[GRAPHSWITCHING_BLOCK_COUNT] = {0};
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
                graphswitching_word_t allowed[GRAPHSWITCHING_BLOCK_COUNT] = {0};
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
        const struct graphswitching_search_context *context, graphswitching_word_t mask[], int vertex,
        int adjacent)
{
        int row_offset = vertex * GRAPHSWITCHING_BLOCK_COUNT;
        int word_index;

        for (word_index = 0; word_index < GRAPHSWITCHING_BLOCK_COUNT; ++word_index) {
                graphswitching_word_t neighbours =
                        context->adjacency[row_offset + word_index];

                mask[word_index] &=
                        adjacent ? neighbours : ~neighbours;
        }
}

static void add_adjacency_pattern(
        const struct graphswitching_search_context *context, graphswitching_word_t destination[],
        int first_partition_index, int count, unsigned int pattern)
{
        graphswitching_word_t matching[GRAPHSWITCHING_BLOCK_COUNT];
        int index;
        int word_index;

        for (word_index = 0; word_index < GRAPHSWITCHING_BLOCK_COUNT; ++word_index) {
                matching[word_index] = UINT32_MAX;
        }
        for (index = 0; index < count; ++index) {
                intersect_with_neighbourhood(
                        context, matching,
                        context->partition[
                                first_partition_index + index].vertex,
                        (int)((pattern >> index) & 1U));
        }
        for (word_index = 0; word_index < GRAPHSWITCHING_BLOCK_COUNT; ++word_index) {
                destination[word_index] |= matching[word_index];
        }
}

static int next_mask_vertex(const graphswitching_word_t mask[], int start, int limit)
{
        int vertex;

        for (vertex = start; vertex < limit; ++vertex) {
                if (((mask[vertex >> GRAPHSWITCHING_BLOCK_SHIFT] >>
                      (vertex & GRAPHSWITCHING_BLOCK_MASK)) &
                     UINT32_C(1)) != 0) {
                        return vertex;
                }
        }

        return limit;
}

static int test_regular_11(const struct graphswitching_search_context *context, int current)
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

static int test_regular_22(const struct graphswitching_search_context *context, int current)
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

static int test_regular_12(const struct graphswitching_search_context *context, int current)
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

static int test_regular_21(const struct graphswitching_search_context *context, int current)
{
        int part_size = context->part_size;
        int selected = current - part_size;

        return context->c1_cross_degrees[
                       context->partition[selected + part_size].vertex] ==
               context->c1_cross_degrees[
                       context->partition[part_size].vertex];
}

static int partition_blocks_can_extend(
        const struct graphswitching_search_context *context, int selected_count)
{
        int part_size = context->part_size;
        int first_selected =
                selected_count < part_size ? selected_count : part_size;
        int second_selected =
                selected_count > part_size
                        ? selected_count - part_size
                        : 0;
        int first_remaining = part_size - first_selected;
        int second_remaining = part_size - second_selected;
        int remaining = 2 * part_size - selected_count;
        int forced = 0;
        int vertex;

        for (vertex = 0; vertex < context->vertex_count; ++vertex) {
                int first_observed = 0;
                int second_observed = 0;
                int interval_low;
                int interval_high;
                int index;
                int possible;

                if (context->in_partition[vertex]) {
                        continue;
                }
                if (first_selected == part_size) {
                        first_observed =
                                context->c1_cross_degrees[vertex];
                } else {
                        for (index = 0; index < first_selected; ++index) {
                                first_observed += graphswitching_adjacency_bit(
                                        context, vertex,
                                        context->partition[index].vertex);
                        }
                }
                for (index = 0; index < second_selected; ++index) {
                        second_observed += graphswitching_adjacency_bit(
                                context, vertex,
                                context->partition[part_size + index].vertex);
                }

                interval_low =
                        first_observed > second_observed
                                ? first_observed
                                : second_observed;
                interval_high =
                        first_observed + first_remaining <
                                        second_observed + second_remaining
                                ? first_observed + first_remaining
                                : second_observed + second_remaining;
                possible = interval_low <= interval_high ||
                           (first_observed == first_selected &&
                            second_observed == 0) ||
                           (first_observed == 0 &&
                            second_observed == second_selected);
                if (!possible && ++forced > remaining) {
                        return 0;
                }
        }
        return 1;
}

static enum graphswitching_result apply_partition_global(
        struct graphswitching_search_context *context)
{
        int neighbour_counts[2 * GRAPHSWITCHING_MAX_VERTICES];
        int any_switch = 0;
        int part_size = context->part_size;
        int vertex;

        for (vertex = 0; vertex < context->vertex_count; ++vertex) {
                int index;

                if (context->in_partition[vertex]) {
                        continue;
                }

                neighbour_counts[2 * vertex] =
                        context->c1_cross_degrees[vertex];
                neighbour_counts[2 * vertex + 1] = 0;
                for (index = 0; index < part_size; ++index) {
                        neighbour_counts[2 * vertex + 1] += graphswitching_adjacency_bit(
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
                        graphswitching_write_graph(context);

                apply_wqh(context, neighbour_counts);
                return result;
        }
}
