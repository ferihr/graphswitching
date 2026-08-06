/* Adjacency-matrix input and matrix/graph6 output. */

#include "graphswitching_internal.h"

#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <stdlib.h>

static int parse_order_header(const char *line, int *vertex_count);
static int line_is_blank(const char *line);
static enum graphswitching_result collect_matrix_bits(
        FILE *input, unsigned char bits[], size_t *bit_count,
        int *declared_vertices);
static enum graphswitching_result write_adjacency_matrix(
        const struct graphswitching_search_context *context);
static enum graphswitching_result write_graph6(
        const struct graphswitching_search_context *context);

enum graphswitching_result graphswitching_read_adjacency_matrix(
        struct graphswitching_search_context *context,
        FILE *input, int requested_vertices)
{
        unsigned char bits[GRAPHSWITCHING_MATRIX_BIT_CAPACITY];
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
                word_index = column >> GRAPHSWITCHING_BLOCK_SHIFT;
                bit_index = column & GRAPHSWITCHING_BLOCK_MASK;
                context->adjacency[row * GRAPHSWITCHING_BLOCK_COUNT + word_index] |=
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
                                if (*bit_count >=
                                    GRAPHSWITCHING_MATRIX_BIT_CAPACITY) {
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
                        if (*bit_count >=
                            GRAPHSWITCHING_MATRIX_BIT_CAPACITY) {
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

enum graphswitching_result graphswitching_write_output_header(
        const struct graphswitching_search_context *context)
{
        if (context->output_format == GRAPHSWITCHING_OUTPUT_GRAPH6) {
                return GRAPHSWITCHING_SUCCESS;
        }
        return fprintf(context->output, "n=%d\n", context->vertex_count) < 0
                       ? GRAPHSWITCHING_OUTPUT_ERROR
                       : GRAPHSWITCHING_SUCCESS;
}

enum graphswitching_result graphswitching_write_graph(
        const struct graphswitching_search_context *context)
{
        if (context->output_format == GRAPHSWITCHING_OUTPUT_GRAPH6) {
                return write_graph6(context);
        }
        return write_adjacency_matrix(context);
}

static enum graphswitching_result write_adjacency_matrix(
        const struct graphswitching_search_context *context)
{
        char row[GRAPHSWITCHING_MAX_VERTICES + 1];
        int row_index;

        for (row_index = 0; row_index < context->vertex_count; ++row_index) {
                int column_index;

                for (column_index = 0;
                     column_index < context->vertex_count;
                     ++column_index) {
                        row[column_index] =
                                graphswitching_adjacency_bit(
                                        context, row_index, column_index)
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

static enum graphswitching_result write_graph6(
        const struct graphswitching_search_context *context)
{
        unsigned char byte = 0;
        int bit_count = 0;
        int column;
        int vertex_count = context->vertex_count;

        if (vertex_count <= 62) {
                if (fputc(vertex_count + 63, context->output) == EOF) {
                        return GRAPHSWITCHING_OUTPUT_ERROR;
                }
        } else {
                int shift;

                if (fputc('~', context->output) == EOF) {
                        return GRAPHSWITCHING_OUTPUT_ERROR;
                }
                for (shift = 12; shift >= 0; shift -= 6) {
                        if (fputc(((vertex_count >> shift) & 0x3f) + 63,
                                  context->output) == EOF) {
                                return GRAPHSWITCHING_OUTPUT_ERROR;
                        }
                }
        }

        for (column = 1; column < vertex_count; ++column) {
                int row;

                for (row = 0; row < column; ++row) {
                        byte = (unsigned char)(
                                (byte << 1) |
                                graphswitching_adjacency_bit(
                                        context, row, column));
                        if (++bit_count == 6) {
                                if (fputc(byte + 63, context->output) == EOF) {
                                        return GRAPHSWITCHING_OUTPUT_ERROR;
                                }
                                byte = 0;
                                bit_count = 0;
                        }
                }
        }
        if (bit_count != 0 &&
            fputc((byte << (6 - bit_count)) + 63,
                  context->output) == EOF) {
                return GRAPHSWITCHING_OUTPUT_ERROR;
        }
        return fputc('\n', context->output) == EOF
                       ? GRAPHSWITCHING_OUTPUT_ERROR
                       : GRAPHSWITCHING_SUCCESS;
}
