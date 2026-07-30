#!/usr/bin/env python3
"""Run reproducible graphswitching benchmarks through the nauty pipeline."""

from __future__ import annotations

import argparse
import csv
import os
import resource
import shlex
import shutil
import statistics
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence


ROOT = Path(__file__).resolve().parent.parent
MANIFEST = ROOT / "benchmarks" / "cases.tsv"


@dataclass(frozen=True)
class BenchmarkCase:
    name: str
    suite: str
    program: str
    arguments: tuple[str, ...]
    input_path: str
    expected_unique: int
    focus: str


@dataclass(frozen=True)
class Timing:
    wall: float
    user: float
    system: float
    unique: int


def positive_integer(text: str) -> int:
    try:
        value = int(text)
    except ValueError as error:
        raise argparse.ArgumentTypeError("must be a positive integer") from error
    if value < 1:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return value


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Benchmark graph switching, nauty conversion, canonical labeling, "
            "and duplicate removal."
        )
    )
    parser.add_argument(
        "cases",
        nargs="*",
        metavar="CASE",
        help="named cases to run; when supplied, --suite is ignored",
    )
    parser.add_argument(
        "--suite",
        choices=("quick", "full", "extrafull", "symmetry"),
        default="quick",
        help=(
            "quick runs Sp(6,2); full adds the extended WQH cases; "
            "extrafull also runs the very slow legacy Sp(4,4) case; "
            "symmetry runs the nauty-enabled orbit-representative cases"
        ),
    )
    parser.add_argument(
        "--runs",
        type=positive_integer,
        default=1,
        help="number of repetitions per case; medians are reported (default: 1)",
    )
    parser.add_argument(
        "--list",
        action="store_true",
        help="list the available cases and exit without requiring nauty",
    )
    parser.add_argument(
        "--nauty-prefix",
        metavar="PREFIX",
        help="use PREFIXamtog and PREFIXlabelg instead of auto-detection",
    )
    return parser.parse_args()


def read_manifest() -> list[BenchmarkCase]:
    cases = []
    with MANIFEST.open(encoding="utf-8", newline="") as manifest:
        rows = csv.reader(manifest, delimiter="\t")
        for line_number, row in enumerate(rows, start=1):
            if not row or row[0].startswith("#"):
                continue
            if len(row) != 7:
                raise ValueError(
                    f"{MANIFEST}:{line_number}: expected 7 tab-separated fields"
                )
            name, suite, program, arguments, input_path, expected, focus = row
            if suite not in ("quick", "full", "extrafull", "symmetry"):
                raise ValueError(
                    f"{MANIFEST}:{line_number}: suite must be quick, full, "
                    "extrafull, or symmetry"
                )
            cases.append(
                BenchmarkCase(
                    name=name,
                    suite=suite,
                    program=program,
                    arguments=tuple(shlex.split(arguments)),
                    input_path=input_path,
                    expected_unique=int(expected),
                    focus=focus,
                )
            )

    names = [case.name for case in cases]
    if len(names) != len(set(names)):
        raise ValueError(f"{MANIFEST}: benchmark names must be unique")
    return cases


def list_cases(cases: Sequence[BenchmarkCase]) -> None:
    widths = {
        "name": max(len("case"), *(len(case.name) for case in cases)),
        "suite": max(len("suite"), *(len(case.suite) for case in cases)),
        "program": max(len("program"), *(len(case.program) for case in cases)),
    }
    print(
        f"{'case':<{widths['name']}}  "
        f"{'suite':<{widths['suite']}}  "
        f"{'program':<{widths['program']}}  focus"
    )
    print(
        f"{'-' * widths['name']}  "
        f"{'-' * widths['suite']}  "
        f"{'-' * widths['program']}  "
        f"{'-' * 5}"
    )
    for case in cases:
        print(
            f"{case.name:<{widths['name']}}  "
            f"{case.suite:<{widths['suite']}}  "
            f"{case.program:<{widths['program']}}  "
            f"{case.focus}"
        )


def select_cases(
    cases: Sequence[BenchmarkCase], requested: Sequence[str], suite: str
) -> list[BenchmarkCase]:
    by_name = {case.name: case for case in cases}
    if requested:
        unknown = [name for name in requested if name not in by_name]
        if unknown:
            choices = ", ".join(sorted(by_name))
            raise ValueError(
                f"unknown benchmark case(s): {', '.join(unknown)}; "
                f"available cases: {choices}"
            )
        return [by_name[name] for name in requested]
    if suite == "quick":
        return [case for case in cases if case.suite == "quick"]
    if suite == "full":
        return [
            case for case in cases if case.suite in ("quick", "full")
        ]
    if suite == "extrafull":
        return [
            case
            for case in cases
            if case.suite in ("quick", "full", "extrafull")
        ]
    return [case for case in cases if case.suite == "symmetry"]


def find_command(command: str) -> str | None:
    path = Path(command)
    if path.parent != Path("."):
        return str(path) if path.is_file() and os.access(path, os.X_OK) else None
    return shutil.which(command)


def parse_override(variable: str) -> list[str] | None:
    value = os.environ.get(variable)
    if not value:
        return None
    command = shlex.split(value)
    if not command:
        raise ValueError(f"{variable} does not name a command")
    if find_command(command[0]) is None:
        raise ValueError(f"{variable} command not found: {command[0]}")
    return command


def resolve_nauty(prefix_option: str | None) -> tuple[list[str], list[str]]:
    amtog_override = parse_override("AMTOG")
    labelg_override = parse_override("LABELG")
    if (amtog_override is None) != (labelg_override is None):
        raise ValueError("set both AMTOG and LABELG, or neither")
    if amtog_override is not None and labelg_override is not None:
        return amtog_override, labelg_override

    prefix = (
        prefix_option
        if prefix_option is not None
        else os.environ.get("NAUTY_PREFIX")
    )
    candidates = [(f"{prefix}amtog", f"{prefix}labelg")] if prefix is not None else [
        ("amtog", "labelg"),
        ("nauty-amtog", "nauty-labelg"),
    ]
    for amtog_name, labelg_name in candidates:
        amtog_path = find_command(amtog_name)
        labelg_path = find_command(labelg_name)
        if amtog_path is not None and labelg_path is not None:
            return [amtog_path], [labelg_path]

    if prefix is not None:
        detail = f"commands {prefix}amtog and {prefix}labelg"
    else:
        detail = "amtog/labelg or nauty-amtog/nauty-labelg"
    raise ValueError(
        f"could not find {detail}; install nauty, set NAUTY_PREFIX, "
        "or set both AMTOG and LABELG"
    )


def validate_case_files(cases: Sequence[BenchmarkCase]) -> None:
    missing = []
    for case in cases:
        program = ROOT / case.program
        input_path = ROOT / case.input_path
        if not program.is_file() or not os.access(program, os.X_OK):
            missing.append(f"{case.program} (run make)")
        if not input_path.is_file():
            missing.append(case.input_path)
    if missing:
        details = ", ".join(sorted(set(missing)))
        raise ValueError(f"missing benchmark input(s): {details}")


def run_pipeline(
    case: BenchmarkCase, amtog: Sequence[str], labelg: Sequence[str]
) -> Timing:
    program = ROOT / case.program
    input_path = ROOT / case.input_path
    commands = [
        [str(program), *case.arguments],
        [*amtog, "-q"],
        [*labelg, "-gqt"],
        ["sort", "-u"],
        ["wc", "-l"],
    ]
    sort_environment = os.environ.copy()
    sort_environment["LC_ALL"] = "C"
    processes: list[subprocess.Popen[bytes]] = []
    usage_before = resource.getrusage(resource.RUSAGE_CHILDREN)
    started = time.perf_counter()

    with input_path.open("rb") as graph_input:
        previous_stdout = None
        for index, command in enumerate(commands):
            standard_input = graph_input if index == 0 else previous_stdout
            environment = sort_environment if command[0] == "sort" else None
            process = subprocess.Popen(
                command,
                stdin=standard_input,
                stdout=subprocess.PIPE,
                env=environment,
            )
            processes.append(process)
            if previous_stdout is not None:
                previous_stdout.close()
            previous_stdout = process.stdout

        output, _ = processes[-1].communicate()
        statuses = [process.wait() for process in processes[:-1]]
        statuses.append(processes[-1].returncode)

    wall = time.perf_counter() - started
    usage_after = resource.getrusage(resource.RUSAGE_CHILDREN)
    failures = [
        (shlex.join(command), status)
        for command, status in zip(commands, statuses)
        if status != 0
    ]
    if failures:
        details = ", ".join(
            f"{command} exited {status}" for command, status in failures
        )
        raise RuntimeError(f"{case.name}: pipeline failed: {details}")

    try:
        unique = int(output.strip())
    except ValueError as error:
        raise RuntimeError(
            f"{case.name}: wc returned an invalid count: {output!r}"
        ) from error

    return Timing(
        wall=wall,
        user=usage_after.ru_utime - usage_before.ru_utime,
        system=usage_after.ru_stime - usage_before.ru_stime,
        unique=unique,
    )


def median_timing(timings: Sequence[Timing]) -> Timing:
    unique_counts = {timing.unique for timing in timings}
    if len(unique_counts) != 1:
        raise RuntimeError(f"inconsistent unique counts across runs: {unique_counts}")
    return Timing(
        wall=statistics.median(timing.wall for timing in timings),
        user=statistics.median(timing.user for timing in timings),
        system=statistics.median(timing.system for timing in timings),
        unique=timings[0].unique,
    )


def print_results(
    results: Sequence[tuple[BenchmarkCase, Timing]], runs: int
) -> None:
    case_width = max(len("case"), *(len(case.name) for case, _ in results))
    program_width = max(
        len("program"), *(len(case.program) for case, _ in results)
    )
    heading = "median " if runs > 1 else ""
    print()
    print(
        f"{'case':<{case_width}}  "
        f"{'program':<{program_width}}  "
        f"{'unique':>6}  {'expected':>8}  "
        f"{heading + 'real':>11}  "
        f"{heading + 'user':>11}  "
        f"{heading + 'sys':>11}"
    )
    print(
        f"{'-' * case_width}  "
        f"{'-' * program_width}  "
        f"{'-' * 6}  {'-' * 8}  "
        f"{'-' * 11}  {'-' * 11}  {'-' * 11}"
    )
    for case, timing in results:
        print(
            f"{case.name:<{case_width}}  "
            f"{case.program:<{program_width}}  "
            f"{timing.unique:>6}  {case.expected_unique:>8}  "
            f"{timing.wall:>10.3f}s  "
            f"{timing.user:>10.3f}s  "
            f"{timing.system:>10.3f}s"
        )


def main() -> int:
    arguments = parse_arguments()
    try:
        available = read_manifest()
        if arguments.list:
            list_cases(available)
            return 0

        selected = select_cases(available, arguments.cases, arguments.suite)
        validate_case_files(selected)
        amtog, labelg = resolve_nauty(arguments.nauty_prefix)
    except (OSError, ValueError) as error:
        print(f"benchmark: {error}", file=sys.stderr)
        return 2

    print(f"amtog:  {shlex.join(amtog)}")
    print(f"labelg: {shlex.join(labelg)}")
    print(f"runs:   {arguments.runs}")

    results = []
    failed = False
    for case in selected:
        timings = []
        for run_number in range(1, arguments.runs + 1):
            print(
                f"Running {case.name} ({run_number}/{arguments.runs})...",
                flush=True,
            )
            try:
                timings.append(run_pipeline(case, amtog, labelg))
            except (OSError, RuntimeError) as error:
                print(f"benchmark: {error}", file=sys.stderr)
                return 2
        try:
            timing = median_timing(timings)
        except RuntimeError as error:
            print(f"benchmark: {case.name}: {error}", file=sys.stderr)
            return 2
        if timing.unique != case.expected_unique:
            failed = True
        results.append((case, timing))

    print_results(results, arguments.runs)
    if failed:
        print(
            "\nbenchmark: one or more isomorphism-class counts did not match",
            file=sys.stderr,
        )
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
