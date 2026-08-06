/* Automorphism groups and symmetry-aware search support. */

#include "graphswitching_internal.h"

#include <stdlib.h>
#include <string.h>

#ifdef GRAPHSWITCHING_WITH_NAUTY
#include <nauty.h>
#include <naugroup.h>
#include <schreier.h>

#define AUTO_SYMMETRY_MIN_GROUP_ORDER 5.0

struct permutation_group {
        int vertex_count;
        size_t generator_count;
        size_t generator_capacity;
        int *generators;
        size_t element_count;
        size_t element_capacity;
        uint16_t *elements;
        int allocation_failed;
        double order_mantissa;
        int order_exponent;
};

struct rank_set {
        uint64_t *slots;
        size_t count;
        size_t capacity;
};

struct subset_orbit {
        uint16_t *subsets;
        size_t count;
        size_t capacity;
        int subset_size;
};

struct graphswitching_symmetry_search {
        struct permutation_group graph_group;
        struct permutation_group c1_stabilizer;
        schreier *graph_schreier;
        permnode *graph_generators;
        schreier *c1_schreier;
        permnode *c1_generators;
        struct rank_set seen_gm[GRAPHSWITCHING_GM_SET_SIZE + 1];
        struct rank_set seen_c1[GRAPHSWITCHING_MAX_PART_SIZE + 1];
        struct rank_set seen_c2[GRAPHSWITCHING_MAX_PART_SIZE + 1];
        uint64_t binomial[GRAPHSWITCHING_MAX_VERTICES + 1]
                         [GRAPHSWITCHING_MAX_PART_SIZE + 1];
};

static enum graphswitching_result choose_gm_orbits_recursively(
        struct graphswitching_search_context *context,
        int switching_set[GRAPHSWITCHING_GM_SET_SIZE],
        int selected_count, const int stabilizer_orbits[]);
static enum graphswitching_result prepare_gm_singleton(
        struct graphswitching_search_context *context, int vertex,
        int stabilizer_orbits[], int *is_representative);
static enum graphswitching_result prepare_raw_subset(
        struct graphswitching_search_context *context,
        const int subset[], int subset_size, struct rank_set *seen,
        int stabilizer_orbits[], int *is_representative);
static enum graphswitching_result compute_automorphism_group(
        const struct graphswitching_search_context *context,
        const struct graphswitching_partition_vertex fixed_subset[],
        int fixed_count, struct permutation_group *group);
static void destroy_permutation_group(struct permutation_group *group);
static int permutation_group_order_exceeds(
        const struct permutation_group *group, double limit);
static int permutation_group_order_at_least(
        const struct permutation_group *group, double limit);
static void permutation_group_orbits(
        const struct permutation_group *group, int orbits[]);
static void initialize_schreier_group(
        const struct permutation_group *group,
        schreier **schreier_group, permnode **generators);
static enum graphswitching_result initialize_automorphism_search(
        struct graphswitching_search_context *context);
static void destroy_automorphism_search(
        struct graphswitching_search_context *context);
static enum graphswitching_result mark_partition_subset(
        struct graphswitching_search_context *context,
        int first_partition_index, struct rank_set *seen,
        int *is_representative);
static enum graphswitching_result prepare_subset_orbit(
        struct graphswitching_search_context *context,
        const struct permutation_group *group,
        const int subset[], int subset_size, struct rank_set *seen,
        int stabilizer_orbits[], int *is_representative);
static uint64_t subset_rank(
        const struct graphswitching_search_context *context,
        const int subset[], int subset_size);
static int rank_set_insert(struct rank_set *set, uint64_t rank);
static int rank_set_grow(struct rank_set *set);
static void rank_set_clear(struct rank_set *set);
static void rank_set_destroy(struct rank_set *set);
static int subset_orbit_append(
        struct subset_orbit *orbit, const int subset[]);
static int subset_orbit_grow(struct subset_orbit *orbit);
static void subset_orbit_destroy(struct subset_orbit *orbit);
static void initialize_vertex_orbits(int orbits[], int vertex_count);
static int vertex_orbit_root(int parents[], int vertex);
static void merge_vertex_orbits(int parents[], int first, int second);
static void finish_vertex_orbits(int parents[], int vertex_count);
static void sort_subset(int subset[], int count);

enum graphswitching_result graphswitching_choose_gm_sets_with_symmetry(
        struct graphswitching_search_context *context)
{
        struct graphswitching_symmetry_search *search = context->symmetry;
        int switching_set[GRAPHSWITCHING_GM_SET_SIZE];
        int root_orbits[GRAPHSWITCHING_MAX_VERTICES];
        int size;

        for (size = 1; size <= GRAPHSWITCHING_GM_SET_SIZE; ++size) {
                rank_set_clear(&search->seen_gm[size]);
        }
        permutation_group_orbits(
                &search->graph_group, root_orbits);
        return choose_gm_orbits_recursively(
                context, switching_set, 0, root_orbits);
}

static enum graphswitching_result choose_gm_orbits_recursively(
        struct graphswitching_search_context *context,
        int switching_set[GRAPHSWITCHING_GM_SET_SIZE],
        int selected_count,
        const int stabilizer_orbits[])
{
        struct graphswitching_symmetry_search *search = context->symmetry;
        unsigned char used_orbit[GRAPHSWITCHING_MAX_VERTICES] = {0};
        graphswitching_word_t completion_mask[GRAPHSWITCHING_BLOCK_COUNT];
        int vertex;

        if (selected_count == GRAPHSWITCHING_GM_SET_SIZE - 1) {
                graphswitching_build_gm_regular_completion_mask(
                        context, switching_set, completion_mask);
        }

        /*
         * Retain subset-orbit representatives under Aut(G), but postpone
         * stabilizer construction when testing the remaining vertices is
         * cheaper. In particular, the regular-completion mask makes a
         * three-set stabilizer unnecessary. Large groups also postpone the
         * two-set stabilizer; small groups filter their cached elements.
         */
        for (vertex = 0; vertex < context->vertex_count; ++vertex) {
                enum graphswitching_result result;
                int child_orbits[GRAPHSWITCHING_MAX_VERTICES];
                int already_selected = 0;
                int child[GRAPHSWITCHING_GM_SET_SIZE];
                int child_count = selected_count + 1;
                int is_representative;
                int orbit = stabilizer_orbits[vertex];
                int index;

                for (index = 0; index < selected_count; ++index) {
                        if (switching_set[index] == vertex) {
                                already_selected = 1;
                                break;
                        }
                }
                if (already_selected || used_orbit[orbit]) {
                        continue;
                }
                used_orbit[orbit] = 1;

                if (child_count == GRAPHSWITCHING_GM_SET_SIZE &&
                    ((completion_mask[vertex >> GRAPHSWITCHING_BLOCK_SHIFT] >>
                      (vertex & GRAPHSWITCHING_BLOCK_MASK)) & UINT32_C(1)) == 0) {
                        continue;
                }

                memcpy(child, switching_set,
                       (size_t)selected_count * sizeof(int));
                child[selected_count] = vertex;
                sort_subset(child, child_count);

                if (child_count == GRAPHSWITCHING_GM_SET_SIZE) {
                        result = graphswitching_apply_gm_set(context, child);
                } else {
                        if (child_count == 1) {
                                result = prepare_gm_singleton(
                                        context, child[0], child_orbits,
                                        &is_representative);
                        } else if (child_count == 3 ||
                                   (child_count == 2 &&
                                    search->graph_group.element_count == 0)) {
                                result = prepare_raw_subset(
                                        context, child, child_count,
                                        &search->seen_gm[child_count],
                                        child_orbits,
                                        &is_representative);
                        } else {
                                result = prepare_subset_orbit(
                                        context, &search->graph_group,
                                        child, child_count,
                                        &search->seen_gm[child_count],
                                        child_orbits,
                                        &is_representative);
                        }
                        if (result != GRAPHSWITCHING_SUCCESS) {
                                return result;
                        }
                        if (!is_representative) {
                                continue;
                        }
                        result = choose_gm_orbits_recursively(
                                context, child, child_count,
                                child_orbits);
                }
                if (result != GRAPHSWITCHING_SUCCESS) {
                        return result;
                }
        }

        return GRAPHSWITCHING_SUCCESS;
}

static enum graphswitching_result prepare_gm_singleton(
        struct graphswitching_search_context *context, int vertex,
        int stabilizer_orbits[], int *is_representative)
{
        struct graphswitching_symmetry_search *search = context->symmetry;
        int fixed[1];
        const int *orbits;
        int inserted;

        inserted = rank_set_insert(
                &search->seen_gm[1],
                subset_rank(context, &vertex, 1));
        if (inserted < 0) {
                return GRAPHSWITCHING_MEMORY_ERROR;
        }
        *is_representative = inserted;
        if (!inserted) {
                return GRAPHSWITCHING_SUCCESS;
        }

        fixed[0] = vertex;
        orbits = getorbits(
                fixed, 1, search->graph_schreier,
                &search->graph_generators,
                context->vertex_count);
        memcpy(stabilizer_orbits, orbits,
               (size_t)context->vertex_count * sizeof(int));
        return GRAPHSWITCHING_SUCCESS;
}

static enum graphswitching_result prepare_raw_subset(
        struct graphswitching_search_context *context,
        const int subset[], int subset_size,
        struct rank_set *seen,
        int stabilizer_orbits[],
        int *is_representative)
{
        int inserted = rank_set_insert(
                seen, subset_rank(context, subset, subset_size));

        if (inserted < 0) {
                return GRAPHSWITCHING_MEMORY_ERROR;
        }
        *is_representative = inserted;
        if (inserted) {
                initialize_vertex_orbits(
                        stabilizer_orbits, context->vertex_count);
        }
        return GRAPHSWITCHING_SUCCESS;
}

static struct permutation_group *generator_target;

static void collect_automorphism_generator(
        int count, int *permutation, int *orbits, int numorbits,
        int stabilizer_vertex, int vertex_count)
{
        struct permutation_group *group = generator_target;
        int *resized;
        size_t new_capacity;

        groupautomproc(
                count, permutation, orbits, numorbits,
                stabilizer_vertex, vertex_count);

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

static void collect_group_element(
        int *permutation, int vertex_count,
        int *abort, void *user_data)
{
        struct permutation_group *group =
                (struct permutation_group *)user_data;
        size_t offset;
        int vertex;

        if (group->element_count == group->element_capacity) {
                size_t new_capacity =
                        group->element_capacity == 0
                                ? 16
                                : 2 * group->element_capacity;
                uint16_t *resized = (uint16_t *)realloc(
                        group->elements,
                        new_capacity * (size_t)vertex_count *
                                sizeof(uint16_t));

                if (resized == NULL) {
                        group->allocation_failed = 1;
                        *abort = 1;
                        return;
                }
                group->elements = resized;
                group->element_capacity = new_capacity;
        }

        offset = group->element_count * (size_t)vertex_count;
        for (vertex = 0; vertex < vertex_count; ++vertex) {
                group->elements[offset + (size_t)vertex] =
                        (uint16_t)permutation[vertex];
        }
        ++group->element_count;
}

static enum graphswitching_result compute_automorphism_group(
        const struct graphswitching_search_context *context,
        const struct graphswitching_partition_vertex fixed_subset[],
        int fixed_count,
        struct permutation_group *group)
{
        DEFAULTOPTIONS_GRAPH(options);
        statsblk statistics;
        graph *nauty_graph;
        grouprec *group_record;
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
                        if (graphswitching_adjacency_bit(context, vertex, neighbour)) {
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
        options.userlevelproc = grouplevelproc;
        options.schreier = TRUE;
        generator_target = group;
        nauty_check(WORDSIZE, set_words, vertex_count, NAUTYVERSIONID);
        densenauty(nauty_graph, lab, ptn, orbits, &options, &statistics,
                   set_words, vertex_count, NULL);
        generator_target = NULL;
        group_record = groupptr(FALSE);

        free(nauty_graph);
        free(lab);
        free(ptn);
        free(orbits);

        if (group->allocation_failed) {
                freegroup(group_record);
                destroy_permutation_group(group);
                return GRAPHSWITCHING_MEMORY_ERROR;
        }
        if (statistics.errstatus != 0) {
                freegroup(group_record);
                destroy_permutation_group(group);
                return GRAPHSWITCHING_INVALID_ARGUMENT;
        }
        group->order_mantissa = statistics.grpsize1;
        group->order_exponent = statistics.grpsize2;
        if (!permutation_group_order_exceeds(group, 4096.0)) {
                makecosetreps(group_record);
                (void)allgroup3(
                        group_record, collect_group_element, group);
        }
        freegroup(group_record);
        if (group->allocation_failed) {
                destroy_permutation_group(group);
                return GRAPHSWITCHING_MEMORY_ERROR;
        }

        return GRAPHSWITCHING_SUCCESS;
}

static void destroy_permutation_group(struct permutation_group *group)
{
        free(group->generators);
        free(group->elements);
        memset(group, 0, sizeof(*group));
}

static int permutation_group_order_exceeds(
        const struct permutation_group *group,
        double limit)
{
        double order = group->order_mantissa;
        int exponent;

        for (exponent = 0;
             exponent < group->order_exponent;
             ++exponent) {
                if (order > limit / 10.0) {
                        return 1;
                }
                order *= 10.0;
        }
        return order > limit;
}

static int permutation_group_order_at_least(
        const struct permutation_group *group,
        double limit)
{
        double order = group->order_mantissa;
        int exponent;

        for (exponent = 0;
             exponent < group->order_exponent;
             ++exponent) {
                if (order >= limit / 10.0) {
                        return 1;
                }
                order *= 10.0;
        }
        return order >= limit;
}

static void permutation_group_orbits(
        const struct permutation_group *group,
        int orbits[])
{
        size_t generator;
        int vertex;

        initialize_vertex_orbits(orbits, group->vertex_count);
        for (generator = 0;
             generator < group->generator_count;
             ++generator) {
                const int *permutation =
                        group->generators +
                        generator * (size_t)group->vertex_count;

                for (vertex = 0;
                     vertex < group->vertex_count;
                     ++vertex) {
                        merge_vertex_orbits(
                                orbits, vertex,
                                permutation[vertex]);
                }
        }
        finish_vertex_orbits(orbits, group->vertex_count);
}

static void initialize_schreier_group(
        const struct permutation_group *group,
        schreier **schreier_group,
        permnode **generators)
{
        size_t generator;

        freeschreier(schreier_group, generators);
        newgroup(schreier_group, generators, group->vertex_count);
        for (generator = 0;
             generator < group->generator_count;
             ++generator) {
                int *permutation =
                        group->generators +
                        generator * (size_t)group->vertex_count;

                (void)addgenerator(
                        schreier_group,
                        generators,
                        permutation, group->vertex_count);
        }
        (void)expandschreier(
                *schreier_group, generators, group->vertex_count);
}

void graphswitching_symmetry_partition_orbits(
        struct graphswitching_search_context *context, int current, int orbits[])
{
        struct graphswitching_symmetry_search *search = context->symmetry;
        schreier *schreier_group;
        permnode **generators;
        int fixed[GRAPHSWITCHING_MAX_PART_SIZE];
        int first;
        int fixed_count;
        const int *stabilizer_orbits;
        int index;

        if (current < context->part_size) {
                first = 0;
                fixed_count = current;
                schreier_group = search->graph_schreier;
                generators = &search->graph_generators;
        } else {
                first = context->part_size;
                fixed_count = current - context->part_size;
                schreier_group = search->c1_schreier;
                generators = &search->c1_generators;
        }

        for (index = 0; index < fixed_count; ++index) {
                fixed[index] =
                        context->partition[first + index].vertex;
        }
        stabilizer_orbits = getorbits(
                fixed, fixed_count, schreier_group, generators,
                context->vertex_count);
        memcpy(orbits, stabilizer_orbits,
               (size_t)context->vertex_count * sizeof(int));
}

static enum graphswitching_result initialize_automorphism_search(
        struct graphswitching_search_context *context)
{
        struct graphswitching_symmetry_search *search = context->symmetry;
        enum graphswitching_result result;
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

        result = compute_automorphism_group(
                context, NULL, 0, &search->graph_group);
        if (result == GRAPHSWITCHING_SUCCESS &&
            search->graph_group.generator_count > 0) {
                initialize_schreier_group(
                        &search->graph_group,
                        &search->graph_schreier,
                        &search->graph_generators);
        }
        return result;
}

static void destroy_automorphism_search(struct graphswitching_search_context *context)
{
        struct graphswitching_symmetry_search *search = context->symmetry;
        int size;

        destroy_permutation_group(&search->graph_group);
        destroy_permutation_group(&search->c1_stabilizer);
        freeschreier(
                &search->graph_schreier, &search->graph_generators);
        freeschreier(
                &search->c1_schreier, &search->c1_generators);
        for (size = 0; size <= GRAPHSWITCHING_GM_SET_SIZE; ++size) {
                rank_set_destroy(&search->seen_gm[size]);
        }
        for (size = 0;
             size <= GRAPHSWITCHING_MAX_PART_SIZE;
             ++size) {
                rank_set_destroy(&search->seen_c1[size]);
                rank_set_destroy(&search->seen_c2[size]);
        }
}

enum graphswitching_result graphswitching_symmetry_prepare_c1(
        struct graphswitching_search_context *context, int *is_representative)
{
        struct graphswitching_symmetry_search *search = context->symmetry;
        enum graphswitching_result result;
        int part_size = context->part_size;

        result = mark_partition_subset(
                context, 0,
                &search->seen_c1[part_size], is_representative);
        if (result != GRAPHSWITCHING_SUCCESS ||
            !*is_representative) {
                return result;
        }

        rank_set_clear(&search->seen_c2[part_size]);
        result = compute_automorphism_group(
                context, context->partition, part_size,
                &search->c1_stabilizer);
        if (result == GRAPHSWITCHING_SUCCESS) {
                initialize_schreier_group(
                        &search->c1_stabilizer,
                        &search->c1_schreier,
                        &search->c1_generators);
        }
        return result;
}

enum graphswitching_result graphswitching_symmetry_prepare_c2(
        struct graphswitching_search_context *context, int *is_representative)
{
        int part_size = context->part_size;

        return mark_partition_subset(
                context, part_size,
                &context->symmetry->seen_c2[part_size],
                is_representative);
}

static enum graphswitching_result mark_partition_subset(
        struct graphswitching_search_context *context,
        int first_partition_index,
        struct rank_set *seen,
        int *is_representative)
{
        int subset[GRAPHSWITCHING_MAX_PART_SIZE];
        int part_size = context->part_size;
        int inserted;
        int index;

        for (index = 0; index < part_size; ++index) {
                subset[index] =
                        context->partition[
                                first_partition_index + index].vertex;
        }
        sort_subset(subset, part_size);
        inserted = rank_set_insert(
                seen, subset_rank(context, subset, part_size));
        if (inserted < 0) {
                return GRAPHSWITCHING_MEMORY_ERROR;
        }
        *is_representative = inserted;
        return GRAPHSWITCHING_SUCCESS;
}

static enum graphswitching_result prepare_subset_orbit(
        struct graphswitching_search_context *context,
        const struct permutation_group *group,
        const int subset[],
        int subset_size,
        struct rank_set *seen,
        int stabilizer_orbits[],
        int *is_representative)
{
        struct subset_orbit orbit = {0};
        int image[GRAPHSWITCHING_MAX_PART_SIZE];
        int inserted;
        int store_representatives = stabilizer_orbits != NULL;
        int vertex_count = group->vertex_count;
        size_t position;
        int vertex;

        if (group->element_count > 0) {
                uint64_t minimum_rank =
                        subset_rank(context, subset, subset_size);
                size_t element;

                for (element = 0;
                     element < group->element_count;
                     ++element) {
                        const uint16_t *permutation =
                                group->elements +
                                element * (size_t)vertex_count;
                        uint64_t image_rank;
                        int index;

                        for (index = 0;
                             index < subset_size;
                             ++index) {
                                image[index] = permutation[subset[index]];
                        }
                        sort_subset(image, subset_size);
                        image_rank = subset_rank(
                                context, image, subset_size);
                        if (image_rank < minimum_rank) {
                                minimum_rank = image_rank;
                        }
                }

                inserted = rank_set_insert(seen, minimum_rank);
                if (inserted < 0) {
                        return GRAPHSWITCHING_MEMORY_ERROR;
                }
                if (!inserted) {
                        *is_representative = 0;
                        return GRAPHSWITCHING_SUCCESS;
                }
                *is_representative = 1;

                if (store_representatives) {
                        initialize_vertex_orbits(
                                stabilizer_orbits, vertex_count);
                        for (element = 0;
                             element < group->element_count;
                             ++element) {
                                const uint16_t *permutation =
                                        group->elements +
                                        element *
                                                (size_t)vertex_count;
                                int stabilizes = 1;
                                int index;

                                for (index = 0;
                                     index < subset_size;
                                     ++index) {
                                        image[index] =
                                                permutation[subset[index]];
                                }
                                sort_subset(image, subset_size);
                                for (index = 0;
                                     index < subset_size;
                                     ++index) {
                                        if (image[index] != subset[index]) {
                                                stabilizes = 0;
                                                break;
                                        }
                                }
                                if (stabilizes) {
                                        for (vertex = 0;
                                             vertex < vertex_count;
                                             ++vertex) {
                                                merge_vertex_orbits(
                                                        stabilizer_orbits,
                                                        vertex,
                                                        permutation[vertex]);
                                        }
                                }
                        }
                        finish_vertex_orbits(
                                stabilizer_orbits, vertex_count);
                }
                return GRAPHSWITCHING_SUCCESS;
        }

        inserted = rank_set_insert(
                seen, subset_rank(context, subset, subset_size));
        if (inserted < 0) {
                return GRAPHSWITCHING_MEMORY_ERROR;
        }
        if (!inserted) {
                *is_representative = 0;
                return GRAPHSWITCHING_SUCCESS;
        }
        *is_representative = 1;

        if (store_representatives) {
                struct graphswitching_partition_vertex fixed_subset[
                        GRAPHSWITCHING_MAX_PART_SIZE] = {{0}};
                struct permutation_group stabilizer = {0};
                enum graphswitching_result result;
                int index;

                for (index = 0; index < subset_size; ++index) {
                        fixed_subset[index].vertex = subset[index];
                }
                result = compute_automorphism_group(
                        context, fixed_subset, subset_size,
                        &stabilizer);
                if (result == GRAPHSWITCHING_SUCCESS) {
                        permutation_group_orbits(
                                &stabilizer, stabilizer_orbits);
                }
                destroy_permutation_group(&stabilizer);
                return result;
        }

        orbit.subset_size = subset_size;
        if (!subset_orbit_append(&orbit, subset)) {
                subset_orbit_destroy(&orbit);
                return GRAPHSWITCHING_MEMORY_ERROR;
        }

        for (position = 0; position < orbit.count; ++position) {
                size_t generator;

                for (generator = 0;
                     generator < group->generator_count;
                     ++generator) {
                        const int *permutation =
                                group->generators +
                                generator * (size_t)vertex_count;
                        uint64_t rank;
                        int index;

                        for (index = 0;
                             index < subset_size;
                             ++index) {
                                image[index] = permutation[
                                        orbit.subsets[
                                                position *
                                                        (size_t)subset_size +
                                                (size_t)index]];
                        }
                        sort_subset(image, subset_size);
                        rank = subset_rank(
                                context, image, subset_size);

                        inserted = rank_set_insert(seen, rank);
                        if (inserted < 0) {
                                subset_orbit_destroy(&orbit);
                                return GRAPHSWITCHING_MEMORY_ERROR;
                        }
                        if (inserted &&
                            !subset_orbit_append(&orbit, image)) {
                                subset_orbit_destroy(&orbit);
                                return GRAPHSWITCHING_MEMORY_ERROR;
                        }
                }
        }

        subset_orbit_destroy(&orbit);
        return GRAPHSWITCHING_SUCCESS;
}

static uint64_t subset_rank(const struct graphswitching_search_context *context,
                            const int subset[], int subset_size)
{
        uint64_t rank = 0;
        int index;

        for (index = 0; index < subset_size; ++index) {
                rank += context->symmetry
                                ->binomial[subset[index]][index + 1];
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

static int subset_orbit_append(
        struct subset_orbit *orbit,
        const int subset[])
{
        size_t offset;
        int index;

        if (orbit->count == orbit->capacity &&
            !subset_orbit_grow(orbit)) {
                return 0;
        }

        offset = orbit->count * (size_t)orbit->subset_size;
        for (index = 0; index < orbit->subset_size; ++index) {
                orbit->subsets[offset + (size_t)index] =
                        (uint16_t)subset[index];
        }
        ++orbit->count;
        return 1;
}

static int subset_orbit_grow(struct subset_orbit *orbit)
{
        size_t new_capacity =
                orbit->capacity == 0 ? 16 : 2 * orbit->capacity;
        uint16_t *new_subsets;

        new_subsets = (uint16_t *)malloc(
                new_capacity * (size_t)orbit->subset_size *
                sizeof(uint16_t));
        if (new_subsets == NULL) {
                return 0;
        }

        if (orbit->count > 0) {
                memcpy(new_subsets, orbit->subsets,
                       orbit->count *
                               (size_t)orbit->subset_size *
                               sizeof(uint16_t));
        }
        free(orbit->subsets);
        orbit->subsets = new_subsets;
        orbit->capacity = new_capacity;
        return 1;
}

static void subset_orbit_destroy(struct subset_orbit *orbit)
{
        free(orbit->subsets);
        memset(orbit, 0, sizeof(*orbit));
}

static void initialize_vertex_orbits(int orbits[], int vertex_count)
{
        int vertex;

        for (vertex = 0; vertex < vertex_count; ++vertex) {
                orbits[vertex] = vertex;
        }
}

static int vertex_orbit_root(int parents[], int vertex)
{
        int root = vertex;

        while (parents[root] != root) {
                root = parents[root];
        }
        while (parents[vertex] != vertex) {
                int parent = parents[vertex];

                parents[vertex] = root;
                vertex = parent;
        }
        return root;
}

static void merge_vertex_orbits(int parents[], int first, int second)
{
        int first_root = vertex_orbit_root(parents, first);
        int second_root = vertex_orbit_root(parents, second);

        if (first_root != second_root) {
                parents[second_root] = first_root;
        }
}

static void finish_vertex_orbits(int parents[], int vertex_count)
{
        int minimum[GRAPHSWITCHING_MAX_VERTICES];
        int vertex;

        for (vertex = 0; vertex < vertex_count; ++vertex) {
                minimum[vertex] = vertex_count;
        }
        for (vertex = 0; vertex < vertex_count; ++vertex) {
                int root = vertex_orbit_root(parents, vertex);

                if (vertex < minimum[root]) {
                        minimum[root] = vertex;
                }
        }
        for (vertex = 0; vertex < vertex_count; ++vertex) {
                parents[vertex] =
                        minimum[vertex_orbit_root(parents, vertex)];
        }
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

enum graphswitching_result graphswitching_symmetry_initialize(
        struct graphswitching_search_context *context)
{
        enum graphswitching_result result;

        context->symmetry = (struct graphswitching_symmetry_search *)calloc(
                1, sizeof(*context->symmetry));
        if (context->symmetry == NULL) {
                return GRAPHSWITCHING_MEMORY_ERROR;
        }
        result = initialize_automorphism_search(context);
        if (result != GRAPHSWITCHING_SUCCESS) {
                graphswitching_symmetry_destroy(context);
                return result;
        }

        if (context->symmetry->graph_group.generator_count == 0) {
                context->use_symmetry = 0;
        } else if (context->symmetry_mode == GRAPHSWITCHING_SYMMETRY_ON) {
                context->use_symmetry = 1;
        } else {
                context->use_symmetry = permutation_group_order_at_least(
                        &context->symmetry->graph_group,
                        AUTO_SYMMETRY_MIN_GROUP_ORDER);
        }
        return GRAPHSWITCHING_SUCCESS;
}

void graphswitching_symmetry_destroy(
        struct graphswitching_search_context *context)
{
        if (context->symmetry == NULL) {
                return;
        }
        destroy_automorphism_search(context);
        free(context->symmetry);
        context->symmetry = NULL;
        context->use_symmetry = 0;
}

enum graphswitching_result graphswitching_symmetry_accept_gm_set(
        struct graphswitching_search_context *context,
        const int switching_set[GRAPHSWITCHING_GM_SET_SIZE],
        int *is_representative)
{
        return prepare_subset_orbit(
                context, &context->symmetry->graph_group,
                switching_set, GRAPHSWITCHING_GM_SET_SIZE,
                &context->symmetry->seen_gm[GRAPHSWITCHING_GM_SET_SIZE],
                NULL, is_representative);
}

void graphswitching_symmetry_fixed_orbits(
        struct graphswitching_search_context *context,
        int current, int orbits[])
{
        const int *stabilizer_orbits = getorbits(
                context->method_vertices, current,
                context->symmetry->graph_schreier,
                &context->symmetry->graph_generators,
                context->vertex_count);

        memcpy(orbits, stabilizer_orbits,
               (size_t)context->vertex_count * sizeof(int));
}

#else

enum graphswitching_result graphswitching_symmetry_initialize(
        struct graphswitching_search_context *context)
{
        context->use_symmetry = 0;
        return GRAPHSWITCHING_FEATURE_UNAVAILABLE;
}

void graphswitching_symmetry_destroy(
        struct graphswitching_search_context *context)
{
        context->use_symmetry = 0;
}

enum graphswitching_result graphswitching_choose_gm_sets_with_symmetry(
        struct graphswitching_search_context *context)
{
        (void)context;
        return GRAPHSWITCHING_FEATURE_UNAVAILABLE;
}

enum graphswitching_result graphswitching_symmetry_accept_gm_set(
        struct graphswitching_search_context *context,
        const int switching_set[GRAPHSWITCHING_GM_SET_SIZE],
        int *is_representative)
{
        (void)context;
        (void)switching_set;
        *is_representative = 0;
        return GRAPHSWITCHING_FEATURE_UNAVAILABLE;
}

enum graphswitching_result graphswitching_symmetry_prepare_c1(
        struct graphswitching_search_context *context,
        int *is_representative)
{
        (void)context;
        *is_representative = 0;
        return GRAPHSWITCHING_FEATURE_UNAVAILABLE;
}

enum graphswitching_result graphswitching_symmetry_prepare_c2(
        struct graphswitching_search_context *context,
        int *is_representative)
{
        (void)context;
        *is_representative = 0;
        return GRAPHSWITCHING_FEATURE_UNAVAILABLE;
}

void graphswitching_symmetry_partition_orbits(
        struct graphswitching_search_context *context,
        int current, int orbits[])
{
        (void)context;
        (void)current;
        (void)orbits;
}

void graphswitching_symmetry_fixed_orbits(
        struct graphswitching_search_context *context,
        int current, int orbits[])
{
        (void)context;
        (void)current;
        (void)orbits;
}

#endif
