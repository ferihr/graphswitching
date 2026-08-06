#!/usr/bin/env python3
"""Benchmark conservative GCC flag profiles and recommend the fastest one."""

from __future__ import annotations

import argparse
import os
import shlex
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path

import run as benchmark


ROOT = Path(__file__).resolve().parent.parent
PROGRAM = "graphswitching"
TUNING_CASES = (
    "sp6-graphswitching-gm",
    "sp6-wqh-p2",
    "sp6-gm3-aut1-gm-sym",
    "bil223-wqh-p3",
    "sp4-4-wqh-sym-p4",
    "witness-gm-p3",
    "witness-gm-p4-sym",
    "witness-ah-p3",
    "witness-ah-p5-sym",
    "witness-gm2-p2",
    "witness-is5-p3-sym",
    "witness-is3-p4",
    "witness-is5-p4-sym",
    "witness-fano",
)


@dataclass(frozen=True)
class Profile:
    name: str
    flags: tuple[str, ...]

    @property
    def cflags(self) -> str:
        return shlex.join(self.flags)


PROFILES = (
    Profile("o2", ("-O2",)),
    Profile("o3", ("-O3",)),
    Profile("default", ("-O3", "-funroll-loops")),
    Profile("native", ("-O3", "-march=native", "-mtune=native")),
    Profile(
        "native-unroll",
        ("-O3", "-march=native", "-mtune=native", "-funroll-loops"),
    ),
    Profile(
        "native-lto",
        (
            "-O3", "-march=native", "-mtune=native", "-funroll-loops",
            "-flto",
        ),
    ),
)


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "cases", nargs="*", metavar="CASE",
        help=(
            "named graphswitching benchmark cases; when supplied, --suite "
            "is ignored"
        ),
    )
    parser.add_argument(
        "--suite",
        choices=("tuning", "short", "full", "extrafull"),
        default="tuning",
        help="benchmark category to use (default: tuning)",
    )
    parser.add_argument(
        "--runs", type=benchmark.positive_integer, default=3,
        help="runs per case and profile (default: 3)",
    )
    parser.add_argument(
        "--compiler", default=os.environ.get("GCC", "gcc"),
        help="GCC command to test (default: $GCC or gcc)",
    )
    parser.add_argument(
        "--nauty-prefix", metavar="PREFIX",
        help="pass PREFIX to the benchmark's nauty command lookup",
    )
    return parser.parse_args()


def command_from_text(text: str, description: str) -> list[str]:
    command = shlex.split(text)
    if not command or shutil.which(command[0]) is None:
        raise ValueError(f"{description} command not found: {text}")
    return command


def select_graphswitching_cases(
    arguments: argparse.Namespace,
) -> list[benchmark.BenchmarkCase]:
    manifest = benchmark.read_manifest()
    requested = (
        arguments.cases
        if arguments.cases or arguments.suite != "tuning"
        else TUNING_CASES
    )
    selected = benchmark.select_cases(
        manifest, requested, arguments.suite
    )
    unsupported = [case.name for case in selected if case.program != PROGRAM]
    if (arguments.cases or arguments.suite == "tuning") and unsupported:
        raise ValueError(
            "GCC tuning only supports graphswitching cases; unsupported: "
            + ", ".join(unsupported)
        )
    cases = [case for case in selected if case.program == PROGRAM]
    if not cases:
        raise ValueError("the selection contains no graphswitching cases")
    return cases


def build(
    make: list[str], compiler: list[str], profile: Profile,
) -> None:
    temporary_parent = Path(tempfile.gettempdir()).resolve()
    if temporary_parent == ROOT or ROOT in temporary_parent.parents:
        temporary_parent = ROOT.parent
    with tempfile.TemporaryDirectory(
        prefix="graphswitching-gcc-", dir=temporary_parent
    ) as temp_dir:
        environment = os.environ.copy()
        environment["TMPDIR"] = temp_dir
        subprocess.run(
            [
                *make, "-B", "-s", "--no-print-directory", PROGRAM,
                f"CC={shlex.join(compiler)}", f"CFLAGS={profile.cflags}",
            ],
            cwd=ROOT,
            env=environment,
            check=True,
        )


def main() -> int:
    arguments = parse_arguments()
    try:
        compiler = command_from_text(arguments.compiler, "compiler")
        make = command_from_text(os.environ.get("MAKE", "make"), "make")
        cases = select_graphswitching_cases(arguments)
        amtog, labelg = benchmark.resolve_nauty(arguments.nauty_prefix)
        version = subprocess.run(
            [*compiler, "--version"], check=True, text=True,
            stdout=subprocess.PIPE,
        ).stdout.splitlines()[0]
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        print(f"tune-gcc: {error}", file=sys.stderr)
        return 2

    print(f"compiler: {version}")
    print(f"program:  {PROGRAM}")
    print(f"category: {arguments.suite if not arguments.cases else 'named cases'}")
    print(f"cases:    {len(cases)}")
    print(f"runs:     {arguments.runs}")

    scores: dict[Profile, float] = {}
    try:
        for profile in PROFILES:
            print(f"\nBuilding {profile.name}: {profile.cflags}")
            build(make, compiler, profile)
            benchmark.validate_case_files(cases)

            total = 0.0
            for case in cases:
                print(f"  {case.name}...", end="", flush=True)
                timings = [
                    benchmark.run_pipeline(case, amtog, labelg)
                    for _ in range(arguments.runs)
                ]
                timing = benchmark.median_timing(timings)
                if timing.unique != case.expected_unique:
                    raise RuntimeError(
                        f"{case.name}: got {timing.unique} unique classes; "
                        f"expected {case.expected_unique}"
                    )
                total += timing.wall
                print(f" {timing.wall:.3f}s")
            scores[profile] = total
            print(f"  total median generator real: {total:.3f}s")
    except (OSError, RuntimeError, subprocess.CalledProcessError) as error:
        print(f"\ntune-gcc: {error}", file=sys.stderr)
        return 2

    ranked = sorted(scores.items(), key=lambda item: item[1])
    winner, winning_time = ranked[0]
    default_time = scores[next(p for p in PROFILES if p.name == "default")]
    name_width = max(len("profile"), *(len(profile.name) for profile, _ in ranked))
    print("\nResults")
    print(
        f"{'profile':<{name_width}}  "
        f"{'generator':>11}  {'vs default':>10}"
    )
    print(f"{'-' * name_width}  {'-' * 11}  {'-' * 10}")
    for profile, total in ranked:
        print(
            f"{profile.name:<{name_width}}  {total:>10.3f}s  "
            f"{default_time / total:>9.3f}x"
        )

    print(f"\nRecommended CFLAGS={winner.cflags!r}")
    if "-march=native" in winner.flags:
        print("This build is optimized for this CPU and may not be portable.")
    print(f"Rebuilding {PROGRAM} with the recommended profile...")
    try:
        build(make, compiler, winner)
    except (OSError, subprocess.CalledProcessError) as error:
        print(f"tune-gcc: final rebuild failed: {error}", file=sys.stderr)
        return 2
    print(
        f"Selected {winner.name}: {winning_time:.3f}s total median real "
        f"({default_time / winning_time:.3f}x the default profile)."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
