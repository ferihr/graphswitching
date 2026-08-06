/* Public graphswitching API and search dispatch. */

#include "graphswitching_internal.h"

#include <string.h>

void graphswitching_options_init(struct graphswitching_options *options)
{
        if (options == NULL) {
                return;
        }

        options->method = GRAPHSWITCHING_METHOD_GM;
        options->vertex_count = 0;
        options->part_size = 2;
        options->symmetry_mode = GRAPHSWITCHING_SYMMETRY_OFF;
        options->use_symmetry = 0;
        options->output_format = GRAPHSWITCHING_OUTPUT_MATRIX;
}

enum graphswitching_result graphswitching_generate_with_options(
        FILE *input,
        FILE *output,
        const struct graphswitching_options *options)
{
        struct graphswitching_search_context context;
        enum graphswitching_result result;

        if (input == NULL || output == NULL || options == NULL ||
            options->method < GRAPHSWITCHING_METHOD_GM ||
            options->method > GRAPHSWITCHING_METHOD_AH10 ||
            options->vertex_count < 0 ||
            options->vertex_count > GRAPHSWITCHING_MAX_VERTICES ||
            options->symmetry_mode < GRAPHSWITCHING_SYMMETRY_OFF ||
            options->symmetry_mode > GRAPHSWITCHING_SYMMETRY_AUTO ||
            (options->use_symmetry != 0 && options->use_symmetry != 1) ||
            options->output_format < GRAPHSWITCHING_OUTPUT_MATRIX ||
            options->output_format > GRAPHSWITCHING_OUTPUT_GRAPH6 ||
            (options->method == GRAPHSWITCHING_METHOD_WQH &&
             (options->part_size <= 0 ||
              options->part_size > GRAPHSWITCHING_MAX_PART_SIZE))) {
                return GRAPHSWITCHING_INVALID_ARGUMENT;
        }

        memset(&context, 0, sizeof(context));
        context.method = options->method;
        context.part_size = options->part_size;
        context.output = output;
        context.output_format = options->output_format;
        context.symmetry_mode = options->symmetry_mode;
        if (context.symmetry_mode == GRAPHSWITCHING_SYMMETRY_OFF &&
            options->use_symmetry) {
                context.symmetry_mode = GRAPHSWITCHING_SYMMETRY_ON;
        }
        if (options->method != GRAPHSWITCHING_METHOD_WQH) {
                context.definition = graphswitching_method_definition(
                        options->method);
        }

        result = graphswitching_read_adjacency_matrix(
                &context, input, options->vertex_count);
        if (result != GRAPHSWITCHING_SUCCESS) {
                return result;
        }

        if (options->method == GRAPHSWITCHING_METHOD_WQH &&
            2 * context.part_size > context.vertex_count) {
                return GRAPHSWITCHING_INVALID_ARGUMENT;
        }
        if (options->method != GRAPHSWITCHING_METHOD_GM &&
            options->method != GRAPHSWITCHING_METHOD_WQH &&
            (context.definition == NULL ||
             context.definition->order > context.vertex_count)) {
                return GRAPHSWITCHING_INVALID_ARGUMENT;
        }

        if (context.symmetry_mode != GRAPHSWITCHING_SYMMETRY_OFF) {
                result = graphswitching_symmetry_initialize(&context);
                if (result != GRAPHSWITCHING_SUCCESS) {
                        graphswitching_symmetry_destroy(&context);
                        return result;
                }
        }

        result = graphswitching_write_output_header(&context);
        if (result == GRAPHSWITCHING_SUCCESS) {
                if (options->method != GRAPHSWITCHING_METHOD_GM &&
                    options->method != GRAPHSWITCHING_METHOD_WQH) {
                        result = graphswitching_choose_fixed_method(&context);
                } else if (options->method == GRAPHSWITCHING_METHOD_GM) {
                        result = graphswitching_choose_gm_sets(&context);
                } else {
                        result = graphswitching_choose_partition(&context, 0);
                }
        }

        graphswitching_symmetry_destroy(&context);
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
                return "selected symmetry search requires nauty support "
                       "in this build";
        case GRAPHSWITCHING_MEMORY_ERROR:
                return "not enough memory";
        default:
                return "unknown error";
        }
}
