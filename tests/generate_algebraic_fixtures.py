#!/usr/bin/env python3
"""Generate and validate algebraically defined graph fixtures."""

from __future__ import annotations

import argparse
import gzip
import itertools
import random
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


def generalized_quadrangle_gq2_4() -> Matrix:
    """Return the point graph of the elliptic quadric Q^-(5,2)."""

    def quadratic_form(vector: int) -> int:
        coordinates = [(vector >> index) & 1 for index in range(6)]
        return (
            coordinates[0] * coordinates[1]
            ^ coordinates[2] * coordinates[3]
            ^ coordinates[4]
            ^ coordinates[4] * coordinates[5]
            ^ coordinates[5]
        )

    points = [
        vector
        for vector in range(1, 1 << 6)
        if quadratic_form(vector) == 0
    ]
    return adjacency_matrix(
        points,
        lambda left, right: (
            quadratic_form(left)
            ^ quadratic_form(right)
            ^ quadratic_form(left ^ right)
        )
        == 0,
    )


def cycle_graph(order: int) -> Matrix:
    return adjacency_matrix(
        list(range(order)),
        lambda left, right: (left - right) % order in (1, order - 1),
    )


def random_graph(order: int, probability: float, seed: int) -> Matrix:
    generator = random.Random(seed)
    matrix = [[0] * order for _ in range(order)]
    for left in range(order):
        for right in range(left + 1, order):
            if generator.random() < probability:
                matrix[left][right] = 1
                matrix[right][left] = 1
    return matrix


def brinkmann_graph() -> Matrix:
    """Return Brinkmann's 4-regular graph on 21 vertices.

    The adjacency list is the GPL-compatible construction distributed by
    SageMath's ``graphs.BrinkmannGraph``.
    """

    upper_neighbors = {
        0: [2, 5, 7, 13],
        1: [3, 6, 7, 8],
        2: [4, 8, 9],
        3: [5, 9, 10],
        4: [6, 10, 11],
        5: [11, 12],
        6: [12, 13],
        7: [15, 20],
        8: [14, 16],
        9: [15, 17],
        10: [16, 18],
        11: [17, 19],
        12: [18, 20],
        13: [14, 19],
        14: [17, 18],
        15: [18, 19],
        16: [19, 20],
        17: [20],
    }
    matrix = [[0] * 21 for _ in range(21)]
    for left, neighbors in upper_neighbors.items():
        for right in neighbors:
            matrix[left][right] = 1
            matrix[right][left] = 1
    return matrix


def complete_graph(order: int) -> Matrix:
    return adjacency_matrix(
        list(range(order)), lambda left, right: left != right
    )


def affine_orthogonal_polar_graph_6_2_minus() -> Matrix:
    """Return VO^-(6,2) on the 64 vectors of GF(2)^6."""

    def quadratic_form(vector: int) -> int:
        coordinates = [(vector >> index) & 1 for index in range(6)]
        return (
            coordinates[0] * coordinates[1]
            ^ coordinates[2] * coordinates[3]
            ^ coordinates[4]
            ^ coordinates[4] * coordinates[5]
            ^ coordinates[5]
        )

    return adjacency_matrix(
        list(range(1 << 6)),
        lambda left, right: quadratic_form(left ^ right) == 0,
    )


def kneser_graph(order: int, subset_size: int) -> Matrix:
    vertices = list(itertools.combinations(range(order), subset_size))
    return adjacency_matrix(
        vertices,
        lambda left, right: not set(left).intersection(right),
    )


def grassmann_graph_2_5_2() -> Matrix:
    """Return the binary Grassmann graph J_2(5,2)."""

    subspaces = {
        tuple(sorted((left, right, left ^ right)))
        for left in range(1, 1 << 5)
        for right in range(left + 1, 1 << 5)
        if left != right
    }
    vertices = sorted(subspaces)
    return adjacency_matrix(
        vertices,
        lambda left, right: bool(set(left).intersection(right)),
    )


def hamming_graph(dimension: int, alphabet_size: int) -> Matrix:
    vertices = list(itertools.product(range(alphabet_size), repeat=dimension))
    return adjacency_matrix(
        vertices,
        lambda left, right: sum(a != b for a, b in zip(left, right)) == 1,
    )


def validate_graph(matrix: Matrix) -> tuple[int, tuple[int, ...]]:
    order = len(matrix)
    if any(len(row) != order for row in matrix):
        raise ValueError("matrix is not square")
    if any(matrix[index][index] for index in range(order)):
        raise ValueError("matrix has a loop")
    if any(
        matrix[row][column] != matrix[column][row]
        for row in range(order)
        for column in range(row + 1, order)
    ):
        raise ValueError("matrix is not symmetric")
    if any(entry not in (0, 1) for row in matrix for entry in row):
        raise ValueError("matrix contains a non-binary entry")
    return order, tuple(sorted({sum(row) for row in matrix}))


def strongly_regular_parameters(matrix: Matrix) -> tuple[int, int, int, int]:
    order, degrees = validate_graph(matrix)
    if len(degrees) != 1:
        raise ValueError("graph is not regular")
    degree = degrees[0]
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
        "generalized-quadrangle-gq2-4.matrix": (
            generalized_quadrangle_gq2_4(),
            (27, 10, 1, 5),
        ),
    }


def paper_fixtures() -> dict[str, tuple[Matrix, int, tuple[int, ...]]]:
    """Return the examples from Table 2 of the Simoens--Van Overberghe paper.

    The preprint does not specify the random seed for G(20, 1/2), so this
    fixture fixes seed 0 for repeatability.
    """

    return {
        "paper-cycle20.matrix": (cycle_graph(20), 20, (2,)),
        "paper-random20-p05-seed0.matrix": (
            random_graph(20, 0.5, 0),
            20,
            (4, 6, 7, 8, 10, 11, 12, 13, 14),
        ),
        "paper-brinkmann.matrix": (brinkmann_graph(), 21, (4,)),
        "paper-complete30.matrix": (complete_graph(30), 30, (29,)),
        "paper-affine-orthogonal-6-2-minus.matrix": (
            affine_orthogonal_polar_graph_6_2_minus(),
            64,
            (27,),
        ),
        "paper-kneser9-3.matrix": (kneser_graph(9, 3), 84, (20,)),
        "paper-grassmann2-5-2.matrix": (
            grassmann_graph_2_5_2(),
            155,
            (42,),
        ),
        "paper-hamming4-6.matrix.gz": (hamming_graph(4, 6), 1296, (20,)),
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

    for name, (matrix, expected_order, expected_degrees) in (
        paper_fixtures().items()
    ):
        actual_order, actual_degrees = validate_graph(matrix)
        if actual_order != expected_order or actual_degrees != expected_degrees:
            print(
                f"{name}: expected order/degrees "
                f"{(expected_order, expected_degrees)}, generated "
                f"{(actual_order, actual_degrees)}",
                file=sys.stderr,
            )
            failed = True
            continue

        path = ROOT / name
        generated = matrix_text(matrix)
        if arguments.write:
            if path.suffix == ".gz":
                path.write_bytes(
                    gzip.compress(
                        generated.encode("ascii"),
                        compresslevel=9,
                        mtime=0,
                    )
                )
            else:
                path.write_text(generated, encoding="ascii")
            print(f"wrote {path.relative_to(ROOT.parent)}")
        elif not path.is_file() or (
            gzip.decompress(path.read_bytes()).decode("ascii")
            if path.suffix == ".gz"
            else path.read_text(encoding="ascii")
        ) != generated:
            print(
                f"{path.relative_to(ROOT.parent)} is stale; "
                "run this script with --write",
                file=sys.stderr,
            )
            failed = True
        else:
            print(
                f"Checked {name} with order/degrees "
                f"{(actual_order, actual_degrees)}: ok"
            )
    return int(failed)


if __name__ == "__main__":
    sys.exit(main())
