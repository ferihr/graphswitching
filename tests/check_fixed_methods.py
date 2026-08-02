#!/usr/bin/env python3
"""Validate the fixed Simoens--Van Overberghe switching catalogues."""

from __future__ import annotations

import argparse
import os
import re
import shlex
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence, TypeVar


ROOT = Path(__file__).resolve().parent.parent
SOURCE = ROOT / "src" / "switching_methods.c"
BENCHMARK_MANIFEST = ROOT / "benchmarks" / "cases.tsv"
Choice = TypeVar("Choice")


@dataclass(frozen=True)
class Method:
    cli_name: str
    parameter: int | None
    source_name: str
    order: int
    denominator: int
    subgraph_count: int


METHODS = (
    Method("gm", 2, "gm4", 4, 2, 4),
    Method("gm", 3, "gm6", 6, 3, 8),
    Method("wqh", 3, "wqh6", 6, 3, 16),
    Method("ah", 3, "ah6", 6, 2, 2),
    Method("is5", 3, "is6", 6, 5, 6),
    Method("fano", None, "fano", 7, 2, 4),
    Method("gm", 4, "gm8", 8, 4, 22),
    Method("gm2", 2, "gm44", 8, 2, 50),
    Method("wqh", 4, "wqh8", 8, 4, 130),
    Method("is3", 4, "is8_level3", 8, 3, 18),
    Method("is5", 4, "is8_level5", 8, 5, 191),
    Method("ah", 5, "ah10", 10, 2, 48),
)

WITNESS_CONFIG = {
    "gm6": ("gm6", 45, 2),
    "ah6": ("ah6", 150, 1),
    "is6": ("is6", 70, 1),
    "fano": ("fano", 75, 1),
    "gm8": ("gm8", 32, 1),
    "gm44": ("gm44", 40, 1),
    "is8_level3": ("is8-level3", 40, 1),
    "is8_level5": ("is8-level5", 30, 1),
    "ah10": ("ah10", 36, 1),
}


def method_arguments(method: Method) -> list[str]:
    arguments = ["--method", method.cli_name]
    if method.parameter is not None:
        arguments.extend(("--part-size", str(method.parameter)))
    return arguments


def method_label(method: Method) -> str:
    if method.parameter is None:
        return method.cli_name
    return f"{method.cli_name} p={method.parameter}"


def benchmark_case_stem(method: Method) -> str:
    parameter = "" if method.parameter is None else f"-p{method.parameter}"
    return f"witness-{method.cli_name}{parameter}"


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data-only",
        action="store_true",
        help=(
            "check matrices, catalogues, and ordinary search without "
            "running symmetry modes"
        ),
    )
    parser.add_argument(
        "--write-witnesses",
        action="store_true",
        help="write the deterministic nonisomorphic benchmark witnesses",
    )
    return parser.parse_args()


def source_array(source: str, method: Method, element_type: str) -> list[int]:
    suffix = "numerator" if element_type == "int8_t" else "subgraphs"
    match = re.search(
        rf"static const {element_type} {method.source_name}_{suffix}"
        rf"\[\] = \{{(.*?)\}};",
        source,
        re.DOTALL,
    )
    if match is None:
        raise ValueError(f"missing {method.source_name}_{suffix}")
    if element_type == "int8_t":
        return [int(value) for value in re.findall(r"-?\d+", match.group(1))]
    return [
        int(value, 0)
        for value in re.findall(r"0x[0-9a-fA-F]+", match.group(1))
    ]


def edge(mask: int, order: int, first: int, second: int) -> int:
    if first > second:
        first, second = second, first
    bit = first * (2 * order - first - 1) // 2 + second - first - 1
    return (mask >> bit) & 1


def adjacency_matrix(mask: int, order: int) -> list[list[int]]:
    return [
        [0 if row == column else edge(mask, order, row, column)
         for column in range(order)]
        for row in range(order)
    ]


def transform_matrix(
    adjacency: list[list[int]], numerator: list[int], order: int,
    denominator: int,
) -> list[list[int]]:
    denominator_squared = denominator * denominator
    transformed = []
    for first in range(order):
        row_values = []
        for second in range(order):
            scaled = sum(
                numerator[row * order + first]
                * adjacency[row][column]
                * numerator[column * order + second]
                for row in range(order)
                for column in range(order)
            )
            if scaled not in (0, denominator_squared):
                raise ValueError("switching subgraph does not transform to 0/1")
            row_values.append(int(scaled == denominator_squared))
        transformed.append(row_values)
    return transformed


def valid_blocks(
    numerator: list[int], order: int, denominator: int
) -> list[tuple[int, int]]:
    blocks = []
    for mask in range(1 << order):
        image = 0
        for column in range(order):
            scaled = sum(
                numerator[row * order + column]
                for row in range(order)
                if (mask >> row) & 1
            )
            if scaled == denominator:
                image |= 1 << column
            elif scaled != 0:
                break
        else:
            blocks.append((mask, image))
    return blocks


def check_orthogonal(
    numerator: list[int], order: int, denominator: int
) -> None:
    for first in range(order):
        for second in range(order):
            product = sum(
                numerator[row * order + first]
                * numerator[row * order + second]
                for row in range(order)
            )
            expected = denominator * denominator if first == second else 0
            if product != expected:
                raise ValueError("switching matrix is not orthogonal")


def validate_method(source: str, method: Method) -> tuple[bytes, int]:
    numerator = source_array(source, method, "int8_t")
    subgraphs = source_array(source, method, "uint64_t")
    if len(numerator) != method.order * method.order:
        raise ValueError("switching matrix has the wrong order")
    if len(subgraphs) != method.subgraph_count:
        raise ValueError(
            f"expected {method.subgraph_count} irreducible subgraphs, "
            f"found {len(subgraphs)}"
        )
    check_orthogonal(numerator, method.order, method.denominator)

    chosen_adjacency = None
    for subgraph in subgraphs:
        adjacency = adjacency_matrix(subgraph, method.order)
        transformed = transform_matrix(
            adjacency, numerator, method.order, method.denominator
        )
        if any(transformed[index][index] for index in range(method.order)):
            raise ValueError("switching subgraph produces a loop")
        if any(
            transformed[row][column] != transformed[column][row]
            for row in range(method.order)
            for column in range(method.order)
        ):
            raise ValueError("switching subgraph produces a directed graph")
        if chosen_adjacency is None and transformed != adjacency:
            chosen_adjacency = adjacency

    blocks = valid_blocks(numerator, method.order, method.denominator)
    if not blocks or blocks[0][0] != 0:
        raise ValueError("switching method has no empty block")

    order = method.order
    if chosen_adjacency is None:
        moved = next((pair for pair in blocks if pair[0] != pair[1]), None)
        if moved is None:
            raise ValueError("switching method has no nontrivial action")
        chosen_adjacency = adjacency_matrix(subgraphs[0], method.order)
        order += 1
        for row in chosen_adjacency:
            row.append(0)
        chosen_adjacency.append([0] * order)
        for vertex in range(method.order):
            value = (moved[0] >> vertex) & 1
            chosen_adjacency[vertex][method.order] = value
            chosen_adjacency[method.order][vertex] = value

    matrix = "".join(
        "".join(str(entry) for entry in row) + "\n"
        for row in chosen_adjacency
    )
    return matrix.encode("ascii"), order


def witness_path(method: Method) -> Path:
    name, _, _ = WITNESS_CONFIG[method.source_name]
    return ROOT / "tests" / f"witness-{name}.matrix"


def witness_order(method: Method) -> int:
    _, order, _ = WITNESS_CONFIG[method.source_name]
    return order


def witness_matrix(source: str, method: Method) -> bytes:
    """Build a deterministic graph with a planted admissible embedding."""
    numerator = source_array(source, method, "int8_t")
    subgraphs = source_array(source, method, "uint64_t")
    blocks = valid_blocks(numerator, method.order, method.denominator)
    moved_blocks = [block for block in blocks if block[0] != block[1]]
    _, order, state = WITNESS_CONFIG[method.source_name]

    def choose(values: Sequence[Choice]) -> Choice:
        nonlocal state
        state = (1664525 * state + 1013904223) & 0xffffffff
        return values[state % len(values)]

    matrix = adjacency_matrix(choose(subgraphs), method.order)
    matrix = [
        row + [0] * (order - method.order)
        for row in matrix
    ] + [
        [0] * order
        for _ in range(order - method.order)
    ]
    chosen_blocks = [
        choose(blocks) for _ in range(order - method.order)
    ]
    if not any(block != image for block, image in chosen_blocks):
        chosen_blocks[0] = choose(moved_blocks)

    for outside, (block, _) in enumerate(
        chosen_blocks, start=method.order
    ):
        for vertex in range(method.order):
            adjacent = (block >> vertex) & 1
            matrix[vertex][outside] = adjacent
            matrix[outside][vertex] = adjacent
    for first in range(method.order, order):
        for second in range(first + 1, order):
            adjacent = choose([0, 1])
            matrix[first][second] = adjacent
            matrix[second][first] = adjacent

    return "".join(
        "".join(str(entry) for entry in row) + "\n"
        for row in matrix
    ).encode("ascii")


def validate_witness_benchmarks(methods: Sequence[Method]) -> None:
    rows = {}
    for line_number, line in enumerate(
        BENCHMARK_MANIFEST.read_text(encoding="utf-8").splitlines(), start=1,
    ):
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) != 7:
            raise ValueError(
                f"{BENCHMARK_MANIFEST.relative_to(ROOT)}:{line_number}: "
                "expected seven tab-separated fields"
            )
        rows[fields[0]] = fields

    for method in methods:
        stem = benchmark_case_stem(method)
        order = witness_order(method)
        expected_arguments = [
            *method_arguments(method), "--vertices", str(order),
        ]
        expected_input = witness_path(method).relative_to(ROOT).as_posix()
        counts = []
        for name, arguments in (
            (stem, expected_arguments),
            (f"{stem}-sym", [*expected_arguments, "--sym"]),
        ):
            row = rows.get(name)
            if row is None:
                raise ValueError(f"missing positive benchmark case {name}")
            _, suite, program, argument_text, input_path, count, _ = row
            if suite != "full" or program != "graphswitching":
                raise ValueError(f"{name}: expected a full graphswitching case")
            if shlex.split(argument_text) != arguments:
                raise ValueError(f"{name}: arguments do not match its method")
            if input_path != expected_input:
                raise ValueError(f"{name}: does not use its witness graph")
            counts.append(int(count))
        if counts[0] <= 0 or counts[0] != counts[1]:
            raise ValueError(
                f"{stem}: ordinary and --sym expected counts must agree "
                "and be positive"
            )


def command_from_environment(name: str, defaults: tuple[str, ...]) -> list[str]:
    override = os.environ.get(name)
    if override:
        command = shlex.split(override)
        if not command:
            raise ValueError(f"{name} does not name a command")
        return command
    for candidate in defaults:
        path = shutil.which(candidate)
        if path is not None:
            return [path]
    raise ValueError(f"could not find {' or '.join(defaults)}")


def canonical_classes(
    graphswitching: list[str], amtog: list[str], labelg: list[str],
    method: Method, mode: str | None, matrix: bytes, order: int,
) -> set[bytes]:
    mode_arguments = [] if mode is None else [mode]
    generated = subprocess.run(
        [
            *graphswitching,
            *method_arguments(method),
            "--vertices", str(order),
            *mode_arguments,
        ],
        input=matrix,
        stdout=subprocess.PIPE,
        check=True,
    ).stdout
    graph6 = subprocess.run(
        [*amtog, "-q"], input=generated, stdout=subprocess.PIPE, check=True
    ).stdout
    canonical = subprocess.run(
        [*labelg, "-gqt"],
        input=graph6,
        stdout=subprocess.PIPE,
        check=True,
    ).stdout
    return set(canonical.splitlines())


def canonical_input(
    amtog: list[str], labelg: list[str], matrix: bytes, order: int,
) -> bytes:
    graph6 = subprocess.run(
        [*amtog, "-q"],
        input=f"n={order}\n".encode("ascii") + matrix,
        stdout=subprocess.PIPE,
        check=True,
    ).stdout
    canonical = subprocess.run(
        [*labelg, "-gqt"],
        input=graph6,
        stdout=subprocess.PIPE,
        check=True,
    ).stdout.splitlines()
    if len(canonical) != 1:
        raise ValueError("could not canonically label witness input")
    return canonical[0]


def main() -> int:
    arguments = parse_arguments()
    try:
        source = SOURCE.read_text(encoding="ascii")
        fixtures = {
            method: validate_method(source, method) for method in METHODS
        }
        witnesses = {
            method: witness_matrix(source, method)
            for method in METHODS
            if method.source_name in WITNESS_CONFIG
        }
        if arguments.write_witnesses:
            for method, matrix in witnesses.items():
                witness_path(method).write_bytes(matrix)
                print(f"Wrote {witness_path(method).relative_to(ROOT)}")
        for method, matrix in witnesses.items():
            if witness_path(method).read_bytes() != matrix:
                raise ValueError(
                    f"{witness_path(method).relative_to(ROOT)} is stale; "
                    "run with --write-witnesses"
                )
            print(f"Checked {method_label(method)} benchmark witness")
        validate_witness_benchmarks(tuple(witnesses))
        print("Checked paired positive-output benchmark cases")
        for method in METHODS:
            print(
                f"Checked {method_label(method)}: "
                f"{method.subgraph_count} irreducible subgraphs"
            )
        if arguments.data_only:
            graphswitching = command_from_environment(
                "GRAPHSWITCHING", (str(ROOT / "graphswitching"),)
            )
            for method in METHODS:
                matrix, order = fixtures[method]
                generated = subprocess.run(
                    [
                        *graphswitching,
                        *method_arguments(method),
                        "--vertices", str(order),
                    ],
                    input=matrix,
                    stdout=subprocess.PIPE,
                    check=True,
                ).stdout
                if generated == f"n={order}\n".encode("ascii"):
                    raise ValueError(
                        f"{method_label(method)}: ordinary search found no "
                        "switch"
                    )
                print(f"Checked {method_label(method)} ordinary search")
            for method, matrix in witnesses.items():
                order = witness_order(method)
                generated = subprocess.run(
                    [
                        *graphswitching,
                        *method_arguments(method),
                        "--vertices", str(order),
                    ],
                    input=matrix,
                    stdout=subprocess.PIPE,
                    check=True,
                ).stdout
                if generated == f"n={order}\n".encode("ascii"):
                    raise ValueError(
                        f"{method_label(method)}: benchmark witness found "
                        "no switch"
                    )
                print(f"Checked {method_label(method)} benchmark search")
            return 0

        graphswitching = command_from_environment(
            "GRAPHSWITCHING", (str(ROOT / "graphswitching"),)
        )
        amtog = command_from_environment("AMTOG", ("amtog", "nauty-amtog"))
        labelg = command_from_environment(
            "LABELG", ("labelg", "nauty-labelg")
        )
        for method in METHODS:
            matrix, order = fixtures[method]
            ordinary = canonical_classes(
                graphswitching, amtog, labelg,
                method, None, matrix, order,
            )
            symmetric = canonical_classes(
                graphswitching, amtog, labelg,
                method, "--sym", matrix, order,
            )
            if not symmetric or symmetric != ordinary:
                raise ValueError(
                    f"{method_label(method)}: ordinary and --sym class "
                    "sets differ"
                )
            print(
                f"Checked {method_label(method)} ordinary/--sym: "
                f"{len(symmetric)} canonical class(es)"
            )
        for method, matrix in witnesses.items():
            order = witness_order(method)
            ordinary = canonical_classes(
                graphswitching, amtog, labelg,
                method, None, matrix, order,
            )
            symmetric = canonical_classes(
                graphswitching, amtog, labelg,
                method, "--sym", matrix, order,
            )
            original = canonical_input(
                amtog, labelg, matrix, order
            )
            if not symmetric or symmetric != ordinary:
                raise ValueError(
                    f"{method_label(method)} benchmark witness: ordinary "
                    "and --sym class sets differ"
                )
            if original in symmetric:
                raise ValueError(
                    f"{method_label(method)} benchmark witness: switched "
                    "class is isomorphic to the input"
                )
            print(
                f"Checked {method_label(method)} benchmark witness: "
                f"{len(symmetric)} new canonical class(es)"
            )
    except (OSError, subprocess.CalledProcessError, ValueError) as error:
        print(f"check-fixed-methods: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
