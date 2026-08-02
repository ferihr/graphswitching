/* Copyright (C) 2020--2023, Ferdinand Ihringer
 *
 * Command-line interface for graph switching.
 */

#include "graphswitching.h"

#include <errno.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define OUTPUT_BUFFER_SIZE (1024U * 1024U)

static int parse_positive_int(const char *text, int *value)
{
        char *end;
        long parsed;

        if (text == NULL) {
                return 0;
        }

        errno = 0;
        end = NULL;
        parsed = strtol(text, &end, 10);
        if (errno == ERANGE || end == text || *end != '\0' ||
            parsed <= 0 || parsed > INT_MAX) {
                return 0;
        }

        *value = (int)parsed;
        return 1;
}

struct command_line {
        struct graphswitching_options switching;
        const char *input_path;
        const char *output_path;
        int part_size_was_set;
        enum command_line_method {
                COMMAND_METHOD_GM,
                COMMAND_METHOD_WQH,
                COMMAND_METHOD_AH,
                COMMAND_METHOD_GM2,
                COMMAND_METHOD_IS3,
                COMMAND_METHOD_IS5,
                COMMAND_METHOD_FANO
        } method;
};

static const char *program_name(const char *path)
{
        const char *separator = strrchr(path, '/');

        return separator == NULL ? path : separator + 1;
}

static void print_help(FILE *stream, const char *program)
{
        fprintf(stream, "Usage: %s [OPTION]...\n", program);
        fprintf(stream,
                "Enumerate graphs obtained by applying one switching "
                "operation to one adjacency matrix.\n");
        fprintf(stream,
                "All options are optional; GM switching and standard "
                "input/output are the defaults.\n\n");
        fprintf(stream, "Options:\n");
        fprintf(stream,
                "  -m METHOD, --method=METHOD\n"
                "                         select a switching method "
                "(default: gm)\n");
        fprintf(stream,
                "  -n N, --vertices=N    use N vertices instead of "
                "inferring the order\n");
        fprintf(stream,
                "  -p P, --part-size=P   set the method parameter P "
                "(see the method table below)\n");
        fprintf(stream,
                "      --sym              use input-graph symmetries to "
                "search switching orbit representatives\n"
                "                         (requires nauty support in this "
                "build; changes output multiplicity/order)\n");
        fprintf(stream,
                "  -i FILE, --input=FILE read FILE instead of standard "
                "input; '-' means stdin\n");
        fprintf(stream,
                "  -o FILE, --output=FILE\n"
                "                         write FILE instead of standard "
                "output; '-' means stdout\n");
        fprintf(stream,
                "  -h, --help             display this complete help and "
                "exit\n");
        fprintf(stream,
                "  -V, --version          display version information and "
                "exit\n\n");

        fprintf(stream,
                "Switching methods (METHOD, accepted P, catalogue name):\n");
        fprintf(stream,
                "  gm    P=2 (default),3,4   GM4, GM6, GM8\n"
                "  wqh   P=1,...,8           WQH(P,P,N-2P); default P=2\n"
                "  ah    P=3,5               AH6, AH10\n"
                "  gm2   P=2 (default)       GM4,4\n"
                "  is3   P=4                 IS8(level 3)\n"
                "  is5   P=3,4               IS6, IS8(level 5)\n"
                "  fano  no P                Fano\n\n"
                "  --part-size is required for ah, is3, and is5, and is "
                "invalid for fano.\n"
                "  WQH also requires 2P <= N. Every method supports "
                "ordinary and --sym search.\n"
                "  --sym computes Aut(G) once and searches stabilizer-"
                "orbit representatives.\n\n");

        fprintf(stream, "Input:\n");
        fprintf(stream,
                "  A square 0/1 adjacency matrix, optionally preceded by "
                "\"n=N\".\n");
        fprintf(stream,
                "  Whitespace is ignored. Without --vertices or an n=N "
                "line, N is inferred.\n\n");

        fprintf(stream, "Output:\n");
        fprintf(stream,
                "  An \"n=N\" line followed by every switched adjacency "
                "matrix found.\n"
                "  --sym preserves output isomorphism classes, not raw "
                "multiplicity or order.\n\n");

        fprintf(stream, "Limits:\n");
        fprintf(stream, "  1 <= N <= %d; 1 <= WQH P <= %d.\n\n",
                GRAPHSWITCHING_MAX_VERTICES,
                GRAPHSWITCHING_MAX_PART_SIZE);

        fprintf(stream, "Examples:\n");
        fprintf(stream, "  %s < graph.matrix\n", program);
        fprintf(stream,
                "  %s --method wqh --part-size 4 --input graph.matrix\n\n",
                program);
        fprintf(stream, "  %s --method gm --sym < graph.matrix\n", program);
        fprintf(stream,
                "  %s --method wqh --part-size 4 --sym < graph.matrix\n\n",
                program);

        fprintf(stream, "Exit status:\n");
        fprintf(stream,
                "  0  success\n"
                "  1  invalid graph input, output failure, or switching "
                "failure\n"
                "  2  invalid command-line usage\n");
}

static void print_usage_hint(FILE *stream, const char *program)
{
        fprintf(stream, "Try '%s --help' for more information.\n", program);
}

static int option_value(int argc, char *argv[], int *index,
                        const char *long_name, char short_name,
                        const char **value)
{
        const char *argument = argv[*index];
        size_t long_length = strlen(long_name);

        if (strcmp(argument, long_name) == 0 ||
            (argument[0] == '-' && argument[1] == short_name &&
             argument[2] == '\0')) {
                if (*index + 1 >= argc) {
                        return -1;
                }
                *value = argv[++(*index)];
                return **value == '\0' ? -1 : 1;
        }

        if (strncmp(argument, long_name, long_length) == 0 &&
            argument[long_length] == '=') {
                *value = argument + long_length + 1;
                return **value == '\0' ? -1 : 1;
        }

        if (argument[0] == '-' && argument[1] == short_name &&
            argument[2] != '\0') {
                *value = argument + 2;
                return 1;
        }

        return 0;
}

static int parse_method(const char *text, enum command_line_method *method)
{
        if (strcmp(text, "gm") == 0) {
                *method = COMMAND_METHOD_GM;
                return 1;
        }
        if (strcmp(text, "wqh") == 0) {
                *method = COMMAND_METHOD_WQH;
                return 1;
        }
        if (strcmp(text, "ah") == 0) {
                *method = COMMAND_METHOD_AH;
                return 1;
        }
        if (strcmp(text, "gm2") == 0) {
                *method = COMMAND_METHOD_GM2;
                return 1;
        }
        if (strcmp(text, "is3") == 0) {
                *method = COMMAND_METHOD_IS3;
                return 1;
        }
        if (strcmp(text, "is5") == 0) {
                *method = COMMAND_METHOD_IS5;
                return 1;
        }
        if (strcmp(text, "fano") == 0) {
                *method = COMMAND_METHOD_FANO;
                return 1;
        }

        return 0;
}

static int resolve_method(struct command_line *command, const char *program)
{
        int parameter = command->switching.part_size;

        switch (command->method) {
        case COMMAND_METHOD_GM:
                if (parameter == 2) {
                        command->switching.method = GRAPHSWITCHING_METHOD_GM;
                } else if (parameter == 3) {
                        command->switching.method = GRAPHSWITCHING_METHOD_GM6;
                } else if (parameter == 4) {
                        command->switching.method = GRAPHSWITCHING_METHOD_GM8;
                } else {
                        fprintf(stderr,
                                "%s: --method gm supports --part-size "
                                "2, 3, or 4\n",
                                program);
                        return 0;
                }
                return 1;
        case COMMAND_METHOD_WQH:
                if (parameter > GRAPHSWITCHING_MAX_PART_SIZE) {
                        fprintf(stderr,
                                "%s: --part-size must not exceed %d\n",
                                program, GRAPHSWITCHING_MAX_PART_SIZE);
                        return 0;
                }
                command->switching.method = GRAPHSWITCHING_METHOD_WQH;
                return 1;
        case COMMAND_METHOD_AH:
                if (parameter == 3) {
                        command->switching.method = GRAPHSWITCHING_METHOD_AH6;
                } else if (parameter == 5) {
                        command->switching.method = GRAPHSWITCHING_METHOD_AH10;
                } else {
                        fprintf(stderr,
                                "%s: --method ah supports --part-size "
                                "3 or 5\n",
                                program);
                        return 0;
                }
                return 1;
        case COMMAND_METHOD_GM2:
                if (parameter != 2) {
                        fprintf(stderr,
                                "%s: --method gm2 requires --part-size 2\n",
                                program);
                        return 0;
                }
                command->switching.method = GRAPHSWITCHING_METHOD_GM44;
                return 1;
        case COMMAND_METHOD_IS3:
                if (parameter != 4) {
                        fprintf(stderr,
                                "%s: --method is3 requires --part-size 4\n",
                                program);
                        return 0;
                }
                command->switching.method =
                        GRAPHSWITCHING_METHOD_IS8_LEVEL3;
                return 1;
        case COMMAND_METHOD_IS5:
                if (parameter == 3) {
                        command->switching.method = GRAPHSWITCHING_METHOD_IS6;
                } else if (parameter == 4) {
                        command->switching.method =
                                GRAPHSWITCHING_METHOD_IS8_LEVEL5;
                } else {
                        fprintf(stderr,
                                "%s: --method is5 supports --part-size "
                                "3 or 4\n",
                                program);
                        return 0;
                }
                return 1;
        case COMMAND_METHOD_FANO:
                if (command->part_size_was_set) {
                        fprintf(stderr,
                                "%s: --part-size is not valid with "
                                "--method fano\n",
                                program);
                        return 0;
                }
                command->switching.method = GRAPHSWITCHING_METHOD_FANO;
                return 1;
        }

        return 0;
}

static int parse_command_line(int argc, char *argv[],
                              struct command_line *command)
{
        const char *program = program_name(argv[0]);
        int index;

        graphswitching_options_init(&command->switching);
        command->input_path = "-";
        command->output_path = "-";
        command->part_size_was_set = 0;
        command->method = COMMAND_METHOD_GM;

        for (index = 1; index < argc; ++index) {
                const char *argument = argv[index];
                const char *value = NULL;
                int matched;

                if (strcmp(argument, "-h") == 0 ||
                    strcmp(argument, "--help") == 0) {
                        print_help(stdout, program);
                        return 1;
                }
                if (strcmp(argument, "-V") == 0 ||
                    strcmp(argument, "--version") == 0) {
                        printf("%s %s\n", program, GRAPHSWITCHING_VERSION);
                        return 1;
                }
                if (strcmp(argument, "--sym") == 0) {
                        command->switching.use_symmetry = 1;
                        continue;
                }

                matched = option_value(argc, argv, &index, "--method",
                                       'm', &value);
                if (matched != 0) {
                        if (matched < 0) {
                                fprintf(stderr,
                                        "%s: option '%s' requires an "
                                        "argument\n",
                                        program, argument);
                                print_usage_hint(stderr, program);
                                return -1;
                        }
                        if (!parse_method(value,
                                          &command->method)) {
                                fprintf(stderr,
                                        "%s: unknown switching method '%s'\n",
                                        program, value);
                                print_usage_hint(stderr, program);
                                return -1;
                        }
                        continue;
                }

                matched = option_value(argc, argv, &index, "--vertices",
                                       'n', &value);
                if (matched != 0) {
                        if (matched < 0 ||
                            !parse_positive_int(
                                    value,
                                    &command->switching.vertex_count)) {
                                fprintf(stderr,
                                        "%s: --vertices requires a positive "
                                        "integer\n",
                                        program);
                                return -1;
                        }
                        continue;
                }

                matched = option_value(argc, argv, &index, "--part-size",
                                       'p', &value);
                if (matched != 0) {
                        if (matched < 0 ||
                            !parse_positive_int(
                                    value,
                                    &command->switching.part_size)) {
                                fprintf(stderr,
                                        "%s: --part-size requires a positive "
                                        "integer\n",
                                        program);
                                return -1;
                        }
                        command->part_size_was_set = 1;
                        continue;
                }

                matched = option_value(argc, argv, &index, "--input",
                                       'i', &value);
                if (matched != 0) {
                        if (matched < 0) {
                                fprintf(stderr,
                                        "%s: --input requires a file name\n",
                                        program);
                                return -1;
                        }
                        command->input_path = value;
                        continue;
                }

                matched = option_value(argc, argv, &index, "--output",
                                       'o', &value);
                if (matched != 0) {
                        if (matched < 0) {
                                fprintf(stderr,
                                        "%s: --output requires a file name\n",
                                        program);
                                return -1;
                        }
                        command->output_path = value;
                        continue;
                }

                fprintf(stderr, "%s: unexpected argument '%s'\n",
                        program, argument);
                print_usage_hint(stderr, program);
                return -1;
        }

        if (command->switching.vertex_count >
            GRAPHSWITCHING_MAX_VERTICES) {
                fprintf(stderr, "%s: --vertices must not exceed %d\n",
                        program, GRAPHSWITCHING_MAX_VERTICES);
                return -1;
        }
        if (!resolve_method(command, program)) {
                return -1;
        }
        if (strcmp(command->input_path, "-") != 0 &&
            strcmp(command->input_path, command->output_path) == 0) {
                fprintf(stderr,
                        "%s: input and output must not be the same file\n",
                        program);
                return -1;
        }

        return 0;
}

int main(int argc, char *argv[])
{
        static char output_buffer[OUTPUT_BUFFER_SIZE];
        struct command_line command;
        enum graphswitching_result result;
        const char *program = program_name(argv[0]);
        FILE *input = stdin;
        FILE *output = stdout;
        int parse_result;
        int status = EXIT_SUCCESS;

        parse_result = parse_command_line(argc, argv, &command);
        if (parse_result > 0) {
                return EXIT_SUCCESS;
        }
        if (parse_result < 0) {
                return 2;
        }

        if (strcmp(command.input_path, "-") != 0) {
                input = fopen(command.input_path, "rb");
                if (input == NULL) {
                        fprintf(stderr, "%s: cannot open '%s' for reading\n",
                                program, command.input_path);
                        return EXIT_FAILURE;
                }
        }
        if (strcmp(command.output_path, "-") != 0) {
                output = fopen(command.output_path, "wb");
                if (output == NULL) {
                        fprintf(stderr, "%s: cannot open '%s' for writing\n",
                                program, command.output_path);
                        if (input != stdin) {
                                (void)fclose(input);
                        }
                        return EXIT_FAILURE;
                }
        }

        /*
         * Matrix enumeration can produce a large amount of output. A larger
         * full buffer reduces the number of writes without changing the
         * output format.
         */
        (void)setvbuf(output, output_buffer, _IOFBF,
                      sizeof(output_buffer));

        result = graphswitching_generate_with_options(
                input, output, &command.switching);
        if (result != GRAPHSWITCHING_SUCCESS) {
                fprintf(stderr, "%s: %s\n", program,
                        graphswitching_result_string(result));
                status = EXIT_FAILURE;
        }

        if (status == EXIT_SUCCESS && fflush(output) == EOF) {
                fprintf(stderr, "%s: could not write output\n", program);
                status = EXIT_FAILURE;
        }

        if (input != stdin && fclose(input) == EOF &&
            status == EXIT_SUCCESS) {
                fprintf(stderr, "%s: could not close input\n", program);
                status = EXIT_FAILURE;
        }
        if (output != stdout && fclose(output) == EOF &&
            status == EXIT_SUCCESS) {
                fprintf(stderr, "%s: could not close output\n", program);
                status = EXIT_FAILURE;
        }

        return status;
}
