/* Internal data for the fixed switching methods from Simoens--Van
 * Overberghe. */

#ifndef GRAPHSWITCHING_METHODS_H
#define GRAPHSWITCHING_METHODS_H

#include "graphswitching.h"

#include <stddef.h>
#include <stdint.h>

#define GRAPHSWITCHING_MAX_METHOD_ORDER 10
#define GRAPHSWITCHING_MAX_METHOD_BLOCKS \
        (1U << GRAPHSWITCHING_MAX_METHOD_ORDER)

struct switching_method_definition {
        enum graphswitching_method method;
        const char *name;
        int order;
        int denominator;
        const int8_t *numerator;
        const uint64_t *irreducible_subgraphs;
        size_t irreducible_subgraph_count;
};

const struct switching_method_definition *
graphswitching_method_definition(enum graphswitching_method method);

#endif
