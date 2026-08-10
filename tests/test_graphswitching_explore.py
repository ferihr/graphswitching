#!/usr/bin/env python3
"""Tests for graphswitching-explore."""

from __future__ import annotations

import contextlib
import importlib.util
import io
import os
import subprocess
import sys
import tempfile
import time
import unittest
from decimal import Decimal
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
SOURCE = ROOT / "tools" / "graphswitching_explore.py"


def load_module():
    specification = importlib.util.spec_from_file_location(
        "graphswitching_explore", SOURCE
    )
    assert specification is not None
    assert specification.loader is not None
    module = importlib.util.module_from_spec(specification)
    sys.modules[specification.name] = module
    specification.loader.exec_module(module)
    return module


EXPLORE = load_module()


class Graph6Tests(unittest.TestCase):
    def test_decodes_cycle_four(self) -> None:
        rows = [row.decode("ascii") for row in EXPLORE.decode_graph6("Cl")]
        self.assertEqual(rows, ["0101", "1010", "0101", "1010"])

    def test_matrix_input_has_header(self) -> None:
        self.assertEqual(
            EXPLORE.matrix_input("Cl"),
            b"n=4\n0101\n1010\n0101\n1010\n",
        )

    def test_rejects_wrong_record_length(self) -> None:
        with self.assertRaises(EXPLORE.ExploreError):
            EXPLORE.decode_graph6("C")

    def test_automorphism_filter(self) -> None:
        bounds = EXPLORE.AutomorphismFilter(
            exact=None, minimum=Decimal(4), maximum=Decimal(8)
        )
        self.assertFalse(bounds.accepts(Decimal(2)))
        self.assertTrue(bounds.accepts(Decimal(4)))
        self.assertTrue(bounds.accepts(Decimal(8)))
        self.assertFalse(bounds.accepts(Decimal(16)))

    def test_symmetry_is_forwarded_for_both_methods(self) -> None:
        gm = EXPLORE.Switching("gm", 2, True)
        wqh = EXPLORE.Switching("wqh", 3, True)

        self.assertEqual(
            gm.arguments(),
            ["--method", "gm", "--part-size", "2", "--sym"],
        )
        self.assertEqual(
            wqh.arguments(),
            ["--method", "wqh", "--part-size", "3", "--sym"],
        )

    def test_method_parameters_are_validated(self) -> None:
        with contextlib.redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit):
                EXPLORE.parse_arguments(
                    ["--output-dir", "out", "--method", "ah"]
                )
            with self.assertRaises(SystemExit):
                EXPLORE.parse_arguments(
                    [
                        "--output-dir", "out", "--method", "fano",
                        "--part-size", "2",
                    ]
                )
        parsed = EXPLORE.parse_arguments(
            [
                "--output-dir", "out", "--method", "is5",
                "--part-size", "4",
            ]
        )
        self.assertEqual(parsed.part_size, 4)

    def test_rejects_nonfinite_time_limit(self) -> None:
        with self.assertRaises(EXPLORE.argparse.ArgumentTypeError):
            EXPLORE.positive_number("nan")
        with self.assertRaises(EXPLORE.argparse.ArgumentTypeError):
            EXPLORE.positive_number("inf")


class DriverTests(unittest.TestCase):
    @staticmethod
    def write_program(path: Path, body: str) -> None:
        path.write_text("#!/usr/bin/env python3\n" + body, encoding="utf-8")
        path.chmod(0o755)

    def make_fake_tools(self, directory: Path) -> dict[str, str]:
        graphswitching = directory / "fake-graphswitching"
        labelg = directory / "fake-labelg"
        dreadnaut = directory / "fake-dreadnaut"

        self.write_program(
            graphswitching,
            "import os, sys, time\n"
            "if sys.argv[-2:] != ['--format', 'graph6']:\n"
            "    sys.exit(2)\n"
            "sys.stdin.buffer.read()\n"
            "time.sleep(float(os.environ.get('FAKE_SWITCHING_DELAY', '0')))\n"
            "sys.stdout.write('Ch\\n')\n",
        )
        self.write_program(
            labelg,
            "import sys\n"
            "data = sys.stdin.buffer.read()\n"
            "sys.stdout.buffer.write(data)\n",
        )
        self.write_program(
            dreadnaut,
            "import sys\n"
            "data = sys.stdin.read()\n"
            "for _ in range(data.count('\\nx\\n')):\n"
            "    print('1 orbit; grpsize=2; 1 gen;')\n",
        )
        return {
            "GRAPHSWITCHING": str(graphswitching),
            "LABELG": str(labelg),
            "DREADNAUT": str(dreadnaut),
        }

    def run_driver(
        self,
        base: Path,
        extra_arguments: list[str],
        extra_environment: dict[str, str] | None = None,
    ) -> subprocess.CompletedProcess[str]:
        tools = self.make_fake_tools(base)
        seeds = base / "seeds.g6"
        seeds.write_text("Cl\nCl\n", encoding="ascii")
        output = base / "output"
        environment = os.environ.copy()
        environment.update(tools)
        if extra_environment is not None:
            environment.update(extra_environment)
        return subprocess.run(
            [
                sys.executable,
                str(SOURCE),
                "--output-dir",
                str(output),
                "--jobs",
                "2",
                "--rounds",
                "3",
                *extra_arguments,
                str(seeds),
            ],
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            env=environment,
            check=False,
        )

    def test_round_files_contain_only_new_graphs(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary)
            result = self.run_driver(base, [])
            self.assertEqual(result.returncode, 0, result.stderr)
            output = base / "output"
            self.assertEqual(
                (output / "initial.g6").read_text(encoding="ascii"), "Cl\n"
            )
            self.assertEqual(
                (output / "round-0001.g6").read_text(encoding="ascii"),
                "Ch\n",
            )
            self.assertEqual(
                (output / "round-0002.g6").read_text(encoding="ascii"), ""
            )
            self.assertEqual(
                (output / "all.g6").read_text(encoding="ascii"), "Ch\nCl\n"
            )
            self.assertFalse((output / ".work").exists())

    def test_automorphism_size_filter(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary)
            result = self.run_driver(base, ["--aut-size", "3"])
            self.assertEqual(result.returncode, 0, result.stderr)
            output = base / "output"
            self.assertEqual(
                (output / "round-0001.g6").read_text(encoding="ascii"), ""
            )
            self.assertEqual(
                (output / "all.g6").read_text(encoding="ascii"), "Cl\n"
            )

    def test_graph_limit_stops_after_new_class(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary)
            result = self.run_driver(base, ["--max-graphs", "2"])
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("graph limit reached", result.stderr)
            output = base / "output"
            self.assertEqual(
                (output / "all.g6").read_text(encoding="ascii"), "Ch\nCl\n"
            )
            self.assertFalse((output / "round-0002.g6").exists())

    def test_time_limit_terminates_active_worker(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary)
            started = time.monotonic()
            result = self.run_driver(
                base,
                ["--time-limit", "0.1"],
                {"FAKE_SWITCHING_DELAY": "5"},
            )
            elapsed = time.monotonic() - started
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertLess(elapsed, 2.0)
            self.assertIn("time limit reached", result.stderr)
            output = base / "output"
            self.assertEqual(
                (output / "round-0001.g6").read_text(encoding="ascii"), ""
            )


if __name__ == "__main__":
    unittest.main()
