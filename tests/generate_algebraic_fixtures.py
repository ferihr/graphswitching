#!/usr/bin/env python3
"""Generate and validate algebraically defined graph fixtures."""

from __future__ import annotations

import argparse
import itertools
import sys
from pathlib import Path
from typing import Callable, Sequence, TypeVar


ROOT = Path(__file__).resolve().parent
Matrix = list[list[int]]
Vertex = TypeVar("Vertex")


def gf4_multiply(left: int, right: int) -> int:
    """Multiply polynomial-basis GF(4) elements modulo x^2 + x + 1."""
    product = 0
    for bit in range(2):
        if right & (1 << bit):
            product ^= left << bit
    if product & 0b1000:
        product ^= 0b1110
    if product & 0b0100:
        product ^= 0b0111
    return product


def gf4_inverse(value: int) -> int:
    if value == 0:
        raise ValueError("zero has no multiplicative inverse")
    for candidate in range(1, 4):
        if gf4_multiply(value, candidate) == 1:
            return candidate
    raise AssertionError("invalid GF(4) element")


def gf4_projective_points(dimension: int) -> list[tuple[int, ...]]:
    points = []
    for vector in itertools.product(range(4), repeat=dimension):
        first = next((entry for entry in vector if entry != 0), None)
        if first is None:
            continue
        inverse = gf4_inverse(first)
        normalized = tuple(gf4_multiply(entry, inverse) for entry in vector)
        if normalized == vector:
            points.append(vector)
    return points


def symplectic_product(
    left: Sequence[int],
    right: Sequence[int],
    multiply: Callable[[int, int], int],
) -> int:
    value = 0
    for index in range(0, len(left), 2):
        value ^= multiply(left[index], right[index + 1])
        value ^= multiply(left[index + 1], right[index])
    return value


def adjacency_matrix(
    vertices: Sequence[Vertex],
    adjacent: Callable[[Vertex, Vertex], bool],
) -> Matrix:
    return [
        [
            int(row != column and adjacent(left, right))
            for column, right in enumerate(vertices)
        ]
        for row, left in enumerate(vertices)
    ]


def symplectic_sp6_2() -> Matrix:
    points = list(itertools.product(range(2), repeat=6))[1:]
    return adjacency_matrix(
        points,
        lambda left, right: symplectic_product(
            left, right, lambda a, b: a & b
        )
        == 0,
    )


def symplectic_sp4_4() -> Matrix:
    points = gf4_projective_points(4)
    return adjacency_matrix(
        points,
        lambda left, right: symplectic_product(
            left, right, gf4_multiply
        )
        == 0,
    )


def bilinear_forms_bil2_2_3() -> Matrix:
    matrices = list(itertools.product(range(3), repeat=4))

    def rank_one(
        left: tuple[int, ...], right: tuple[int, ...]
    ) -> bool:
        difference = tuple(
            (a - b) % 3
            for a, b in zip(left, right)
        )
        determinant = (
            difference[0] * difference[3]
            - difference[1] * difference[2]
        ) % 3
        return determinant == 0 and any(difference)

    return adjacency_matrix(matrices, rank_one)


def strongly_regular_parameters(matrix: Matrix) -> tuple[int, int, int, int]:
    order = len(matrix)
    if any(len(row) != order for row in matrix):
        raise ValueError("matrix is not square")
    if any(matrix[index][index] for index in range(order)):
        raise ValueError("matrix has a loop")
    if any(
        matrix[row][column] != matrix[column][row]
        for row in range(order)
        for column in range(order)
    ):
        raise ValueError("matrix is not symmetric")

    degrees = {sum(row) for row in matrix}
    if len(degrees) != 1:
        raise ValueError("graph is not regular")
    degree = degrees.pop()
    adjacent_common = set()
    nonadjacent_common = set()
    for left in range(order):
        for right in range(left + 1, order):
            common = sum(
                matrix[left][vertex] and matrix[right][vertex]
                for vertex in range(order)
            )
            target = (
                adjacent_common
                if matrix[left][right]
                else nonadjacent_common
            )
            target.add(common)
    if len(adjacent_common) != 1 or len(nonadjacent_common) != 1:
        raise ValueError("graph is not strongly regular")
    return (
        order,
        degree,
        adjacent_common.pop(),
        nonadjacent_common.pop(),
    )


def matrix_text(matrix: Matrix) -> str:
    return "".join("".join(map(str, row)) + "\n" for row in matrix)


def fixtures() -> dict[str, tuple[Matrix, tuple[int, int, int, int]]]:
    return {
        "symplectic-sp6-2.matrix": (
            symplectic_sp6_2(),
            (63, 30, 13, 15),
        ),
        "bilinear-forms-bil2-2-3.matrix": (
            bilinear_forms_bil2_2_3(),
            (81, 32, 13, 12),
        ),
        "symplectic-sp4-4.matrix": (
            symplectic_sp4_4(),
            (85, 20, 3, 5),
        ),
    }


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="generate or validate the algebraic graph fixtures"
    )
    parser.add_argument(
        "--write",
        action="store_true",
        help="write the generated matrices; the default only checks them",
    )
    return parser.parse_args()


def main() -> int:
    arguments = parse_arguments()
    failed = False
    for name, (matrix, expected_parameters) in fixtures().items():
        actual_parameters = strongly_regular_parameters(matrix)
        if actual_parameters != expected_parameters:
            print(
                f"{name}: expected parameters {expected_parameters}, "
                f"generated {actual_parameters}",
                file=sys.stderr,
            )
            failed = True
            continue

        path = ROOT / name
        generated = matrix_text(matrix)
        if arguments.write:
            path.write_text(generated, encoding="ascii")
            print(f"wrote {path.relative_to(ROOT.parent)}")
        elif not path.is_file() or path.read_text(encoding="ascii") != generated:
            print(
                f"{path.relative_to(ROOT.parent)} is stale; "
                "run this script with --write",
                file=sys.stderr,
            )
            failed = True
        else:
            print(
                f"Checked {name} as SRG{actual_parameters}: ok"
            )
    return int(failed)


if __name__ == "__main__":
    sys.exit(main())
