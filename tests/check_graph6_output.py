#!/usr/bin/env python3
"""Check direct graph6 output against the default matrix records."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent


def matrix_records(output: bytes) -> list[bytes]:
    lines = output.splitlines()
    if not lines or not lines[0].startswith(b"n="):
        raise ValueError("matrix output has no order header")
    order = int(lines[0][2:])
    rows = [line for line in lines[1:] if line]
    if len(rows) % order:
        raise ValueError("matrix output has an incomplete record")
    return [
        b"\n".join(rows[index:index + order]) + b"\n"
        for index in range(0, len(rows), order)
    ]


def decode_graph6(record: bytes) -> bytes:
    values = [value - 63 for value in record.strip()]
    if not values:
        raise ValueError("empty graph6 record")
    if values[0] <= 62:
        order = values[0]
        payload = values[1:]
    elif len(values) >= 4 and values[0] == 63:
        order = (values[1] << 12) | (values[2] << 6) | values[3]
        payload = values[4:]
    else:
        raise ValueError("unsupported graph6 order header")

    bits = [
        (value >> shift) & 1
        for value in payload
        for shift in range(5, -1, -1)
    ]
    matrix = [[0] * order for _ in range(order)]
    position = 0
    for column in range(1, order):
        for row in range(column):
            matrix[row][column] = bits[position]
            matrix[column][row] = bits[position]
            position += 1
    return b"".join(
        bytes(48 + value for value in row) + b"\n" for row in matrix
    )


def check_case(arguments: list[str], input_path: str) -> None:
    matrix = (ROOT / input_path).read_bytes()
    command = [str(ROOT / "graphswitching"), *arguments]
    ordinary = subprocess.run(
        command, input=matrix, stdout=subprocess.PIPE, check=True
    ).stdout
    graph6 = subprocess.run(
        [*command, "--format", "graph6"],
        input=matrix,
        stdout=subprocess.PIPE,
        check=True,
    ).stdout
    expected = matrix_records(ordinary)
    actual = [decode_graph6(record) for record in graph6.splitlines()]
    if actual != expected:
        raise ValueError(f"{' '.join(arguments)}: graph6 records differ")


def main() -> int:
    try:
        check_case(["--vertices", "10"], "tests/petersen.matrix")
        check_case(
            ["--method", "wqh", "--vertices", "63", "--part-size", "2"],
            "tests/symplectic-sp6-2.matrix",
        )
    except (OSError, subprocess.CalledProcessError, ValueError) as error:
        print(f"check-graph6-output: {error}", file=sys.stderr)
        return 1
    print("Direct graph6 output matches matrix output.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
