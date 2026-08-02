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
        help="named benchmark cases; when supplied, --suite is ignored",
    )
    parser.add_argument(
        "--suite", choices=("short", "full", "extrafull"), default="short",
        help="benchmark suite to use (default: short)",
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


def build(
    make: list[str], compiler: list[str], profile: Profile,
    programs: list[str],
) -> None:
    with tempfile.TemporaryDirectory(
        prefix=".graphswitching-gcc-", dir=ROOT
    ) as temp_dir:
        environment = os.environ.copy()
        environment["TMPDIR"] = temp_dir
        subprocess.run(
            [
                *make, "-B", "-s", "--no-print-directory", *programs,
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
        cases = benchmark.select_cases(
            benchmark.read_manifest(), arguments.cases, arguments.suite
        )
        amtog, labelg = benchmark.resolve_nauty(arguments.nauty_prefix)
        version = subprocess.run(
            [*compiler, "--version"], check=True, text=True,
            stdout=subprocess.PIPE,
        ).stdout.splitlines()[0]
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        print(f"tune-gcc: {error}", file=sys.stderr)
        return 2

    programs = list(dict.fromkeys(case.program for case in cases))
    print(f"compiler: {version}")
    print(f"suite:    {arguments.suite if not arguments.cases else 'named cases'}")
    print(f"cases:    {len(cases)}")
    print(f"runs:     {arguments.runs}")

    scores: dict[Profile, float] = {}
    try:
        for profile in PROFILES:
            print(f"\nBuilding {profile.name}: {profile.cflags}")
            build(make, compiler, profile, programs)
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
            print(f"  total median real: {total:.3f}s")
    except (OSError, RuntimeError, subprocess.CalledProcessError) as error:
        print(f"\ntune-gcc: {error}", file=sys.stderr)
        return 2

    ranked = sorted(scores.items(), key=lambda item: item[1])
    winner, winning_time = ranked[0]
    default_time = scores[next(p for p in PROFILES if p.name == "default")]
    name_width = max(len("profile"), *(len(profile.name) for profile, _ in ranked))
    print("\nResults")
    print(f"{'profile':<{name_width}}  {'median real':>11}  {'vs default':>10}")
    print(f"{'-' * name_width}  {'-' * 11}  {'-' * 10}")
    for profile, total in ranked:
        print(
            f"{profile.name:<{name_width}}  {total:>10.3f}s  "
            f"{default_time / total:>9.3f}x"
        )

    print(f"\nRecommended CFLAGS={winner.cflags!r}")
    if "-march=native" in winner.flags:
        print("This build is optimized for this CPU and may not be portable.")
    print("Rebuilding the benchmark programs with the recommended profile...")
    try:
        build(make, compiler, winner, programs)
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
