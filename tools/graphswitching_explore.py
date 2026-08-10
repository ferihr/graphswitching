#!/usr/bin/env python3
"""Explore the switching closure of graph6 seed graphs."""

from __future__ import annotations

import argparse
import concurrent.futures
import dataclasses
import decimal
import math
import os
import re
import shlex
import shutil
import signal
import subprocess
import sys
import threading
import time
from pathlib import Path
from typing import Iterable, Sequence


VERSION = "@GRAPHSWITCHING_VERSION@"
GRAPH6_HEADER = ">>graph6<<"
MAX_VERTICES = 967
SWITCHING_METHODS = (
    "gm",
    "wqh",
    "ah",
    "gm2",
    "is3",
    "is5",
    "fano",
)
AUTOMORPHISM_SIZE_PATTERN = re.compile(r"grpsize=([^;]+);")
AUTOMORPHISM_BATCH_SIZE = 64


class ExploreError(Exception):
    """A user-facing exploration error."""


@dataclasses.dataclass(frozen=True)
class Toolchain:
    graphswitching: tuple[str, ...]
    labelg: tuple[str, ...]
    dreadnaut: tuple[str, ...] | None
    sort: tuple[str, ...]


@dataclasses.dataclass(frozen=True)
class Switching:
    method: str
    part_size: int
    symmetry: bool

    def arguments(self) -> list[str]:
        arguments = ["--method", self.method]
        if self.method != "fano":
            arguments.extend(("--part-size", str(self.part_size)))
        if self.symmetry:
            arguments.append("--sym")
        return arguments


@dataclasses.dataclass(frozen=True)
class AutomorphismFilter:
    exact: decimal.Decimal | None
    minimum: decimal.Decimal | None
    maximum: decimal.Decimal | None

    def accepts(self, size: decimal.Decimal) -> bool:
        if self.exact is not None and size != self.exact:
            return False
        if self.minimum is not None and size < self.minimum:
            return False
        if self.maximum is not None and size > self.maximum:
            return False
        return True


@dataclasses.dataclass(frozen=True)
class JobResult:
    output_paths: tuple[Path, ...]
    error_path: Path


@dataclasses.dataclass
class RoundResult:
    graphs: set[str]
    processed: int
    timed_out: bool
    graph_limit_reached: bool


def program_name() -> str:
    return Path(sys.argv[0]).name


def positive_integer(text: str) -> int:
    try:
        value = int(text)
    except ValueError as error:
        raise argparse.ArgumentTypeError("must be a positive integer") from error
    if value < 1:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return value


def nonnegative_integer(text: str) -> int:
    try:
        value = int(text)
    except ValueError as error:
        raise argparse.ArgumentTypeError(
            "must be a nonnegative integer"
        ) from error
    if value < 0:
        raise argparse.ArgumentTypeError("must be a nonnegative integer")
    return value


def positive_number(text: str) -> float:
    try:
        value = float(text)
    except ValueError as error:
        raise argparse.ArgumentTypeError("must be a positive number") from error
    if not math.isfinite(value) or value <= 0.0:
        raise argparse.ArgumentTypeError("must be a positive number")
    return value


def positive_decimal_integer(text: str) -> decimal.Decimal:
    if (
        not text.isascii()
        or not text.isdigit()
        or not any(character != "0" for character in text)
    ):
        raise argparse.ArgumentTypeError("must be a positive decimal integer")
    return decimal.Decimal(text)


def parse_arguments(arguments: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Repeatedly apply graphswitching to graph6 graphs, retaining one "
            "canonical representative of each new isomorphism class."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Output files:\n"
            "  initial.g6       canonical unique seeds\n"
            "  round-NNNN.g6    classes first found in that round\n"
            "  all.g6           accumulated classes\n\n"
            "Environment overrides:\n"
            "  GRAPHSWITCHING, LABELG, DREADNAUT, SORT\n"
        ),
    )
    parser.add_argument(
        "inputs",
        nargs="*",
        metavar="FILE",
        help=(
            "graph6 input files; use '-' for standard input "
            "(default: standard input)"
        ),
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        required=True,
        metavar="DIR",
        help="write every generated file below the new or empty directory DIR",
    )
    parser.add_argument(
        "-j",
        "--jobs",
        "--processes",
        type=positive_integer,
        default=max(1, os.cpu_count() or 1),
        metavar="N",
        help=(
            "run up to N switching jobs and apportion N canonical-labelling "
            "lanes among them (default: CPUs)"
        ),
    )
    parser.add_argument(
        "-r",
        "--rounds",
        type=nonnegative_integer,
        default=10,
        metavar="N",
        help="stop after N switching rounds (default: 10)",
    )
    parser.add_argument(
        "--time-limit",
        type=positive_number,
        metavar="SECONDS",
        help="stop after this many seconds of switching work",
    )
    parser.add_argument(
        "--max-graphs",
        type=positive_integer,
        metavar="N",
        help="stop when initial.g6 and the round files contain N graphs",
    )
    parser.add_argument(
        "-m",
        "--method",
        choices=SWITCHING_METHODS,
        default="gm",
        help="switching method passed to graphswitching (default: gm)",
    )
    parser.add_argument(
        "-p",
        "--part-size",
        type=positive_integer,
        default=None,
        metavar="P",
        help="method parameter passed to graphswitching (gm/wqh default: 2)",
    )
    parser.add_argument(
        "--sym",
        action="store_true",
        help="pass --sym to graphswitching",
    )
    automorphisms = parser.add_argument_group(
        "automorphism group size filter"
    )
    automorphisms.add_argument(
        "--aut-size",
        type=positive_decimal_integer,
        metavar="N",
        help="keep a newly generated graph only when |Aut(G)| equals N",
    )
    automorphisms.add_argument(
        "--min-aut-size",
        type=positive_decimal_integer,
        metavar="N",
        help="keep a newly generated graph only when |Aut(G)| is at least N",
    )
    automorphisms.add_argument(
        "--max-aut-size",
        type=positive_decimal_integer,
        metavar="N",
        help="keep a newly generated graph only when |Aut(G)| is at most N",
    )
    programs = parser.add_argument_group("program locations")
    programs.add_argument(
        "--graphswitching",
        metavar="COMMAND",
        help="graphswitching command (default: sibling executable or PATH)",
    )
    programs.add_argument(
        "--nauty-prefix",
        metavar="PREFIX",
        help="use PREFIXlabelg and PREFIXdreadnaut",
    )
    programs.add_argument(
        "--labelg",
        metavar="COMMAND",
        help="labelg command, overriding LABELG and auto-detection",
    )
    programs.add_argument(
        "--dreadnaut",
        metavar="COMMAND",
        help="dreadnaut command used by automorphism filters",
    )
    parser.add_argument(
        "-q",
        "--quiet",
        action="store_true",
        help="suppress per-round progress messages",
    )
    parser.add_argument(
        "-V",
        "--version",
        action="version",
        version=f"%(prog)s {VERSION}",
    )
    parsed = parser.parse_args(arguments)

    supplied_parameter = parsed.part_size
    if supplied_parameter is None:
        if parsed.method in ("gm", "wqh", "gm2", "fano"):
            parsed.part_size = 2
        else:
            parser.error(f"--method {parsed.method} requires --part-size")
    supported_parameters = {
        "gm": (2, 3, 4),
        "ah": (3, 5),
        "gm2": (2,),
        "is3": (4,),
        "is5": (3, 4),
    }
    if parsed.method == "fano" and supplied_parameter is not None:
        parser.error("--part-size is not valid with --method fano")
    if (
        parsed.method in supported_parameters
        and parsed.part_size not in supported_parameters[parsed.method]
    ):
        choices = ", ".join(
            str(value) for value in supported_parameters[parsed.method]
        )
        parser.error(
            f"--method {parsed.method} supports --part-size {choices}"
        )
    if parsed.method == "wqh" and parsed.part_size > 8:
        parser.error("--method wqh supports --part-size at most 8")

    if parsed.aut_size is not None and (
        parsed.min_aut_size is not None or parsed.max_aut_size is not None
    ):
        parser.error("--aut-size cannot be combined with an automorphism bound")
    if (
        parsed.min_aut_size is not None
        and parsed.max_aut_size is not None
        and parsed.min_aut_size > parsed.max_aut_size
    ):
        parser.error("--min-aut-size must not exceed --max-aut-size")
    return parsed


def split_command(text: str, description: str) -> tuple[str, ...]:
    command = tuple(shlex.split(text))
    if not command:
        raise ExploreError(f"{description} does not name a command")
    executable = command[0]
    if os.sep in executable:
        path = Path(executable).expanduser()
        if not path.is_file() or not os.access(path, os.X_OK):
            raise ExploreError(f"{description} is not executable: {executable}")
        command = (str(path), *command[1:])
    elif shutil.which(executable) is None:
        raise ExploreError(f"{description} command not found: {executable}")
    return command


def resolve_command(
    explicit: str | None,
    environment_name: str,
    candidates: Iterable[str],
    description: str,
) -> tuple[str, ...]:
    if explicit is not None:
        return split_command(explicit, description)
    environment_value = os.environ.get(environment_name)
    if environment_value:
        return split_command(environment_value, environment_name)
    for candidate in candidates:
        if os.sep in candidate:
            path = Path(candidate)
            if path.is_file() and os.access(path, os.X_OK):
                return (str(path),)
        else:
            resolved = shutil.which(candidate)
            if resolved is not None:
                return (resolved,)
    names = ", ".join(candidates)
    raise ExploreError(f"could not find {description}; tried {names}")


def resolve_toolchain(options: argparse.Namespace) -> Toolchain:
    invoked_dir = Path(sys.argv[0]).resolve().parent
    source_dir = Path(__file__).resolve().parent
    graphswitching_candidates = (
        str(invoked_dir / "graphswitching"),
        str(source_dir / "graphswitching"),
        str(source_dir.parent / "graphswitching"),
        "graphswitching",
    )
    graphswitching = resolve_command(
        options.graphswitching,
        "GRAPHSWITCHING",
        graphswitching_candidates,
        "graphswitching",
    )

    if options.nauty_prefix is not None:
        labelg_candidates = (f"{options.nauty_prefix}labelg",)
        dreadnaut_candidates = (f"{options.nauty_prefix}dreadnaut",)
    else:
        labelg_candidates = ("labelg", "nauty-labelg")
        dreadnaut_candidates = ("dreadnaut", "nauty-dreadnaut")

    labelg = resolve_command(
        options.labelg, "LABELG", labelg_candidates, "nauty labelg"
    )
    automorphism_filter_requested = any(
        value is not None
        for value in (
            options.aut_size,
            options.min_aut_size,
            options.max_aut_size,
        )
    )
    dreadnaut = None
    if automorphism_filter_requested:
        dreadnaut = resolve_command(
            options.dreadnaut,
            "DREADNAUT",
            dreadnaut_candidates,
            "nauty dreadnaut",
        )
    sort = resolve_command(
        None, "SORT", ("sort",), "sort"
    )
    return Toolchain(
        graphswitching=graphswitching,
        labelg=labelg,
        dreadnaut=dreadnaut,
        sort=sort,
    )


def graph6_order_and_offset(data: bytes) -> tuple[int, int]:
    if not data:
        raise ExploreError("empty graph6 record")
    values = [byte - 63 for byte in data]
    if any(value < 0 or value > 63 for value in values):
        raise ExploreError("graph6 records must contain ASCII bytes 63 through 126")
    if values[0] < 63:
        return values[0], 1
    if len(values) < 4:
        raise ExploreError("truncated graph6 order")
    if values[1] < 63:
        order = (values[1] << 12) | (values[2] << 6) | values[3]
        return order, 4
    if len(values) < 8:
        raise ExploreError("truncated large graph6 order")
    order = 0
    for value in values[2:8]:
        order = (order << 6) | value
    return order, 8


def decode_graph6(record: str) -> list[bytearray]:
    try:
        data = record.encode("ascii")
    except UnicodeEncodeError as error:
        raise ExploreError("graph6 records must be ASCII") from error
    order, offset = graph6_order_and_offset(data)
    if order < 1:
        raise ExploreError("graphswitching requires at least one vertex")
    if order > MAX_VERTICES:
        raise ExploreError(
            f"graphswitching supports at most {MAX_VERTICES} vertices"
        )
    bit_count = order * (order - 1) // 2
    data_count = (bit_count + 5) // 6
    if len(data) != offset + data_count:
        raise ExploreError(
            f"graph6 record for {order} vertices has the wrong length"
        )

    rows = [bytearray(b"0" * order) for _ in range(order)]
    edge_index = 0
    for column in range(1, order):
        for row in range(column):
            value = data[offset + edge_index // 6] - 63
            bit = (value >> (5 - edge_index % 6)) & 1
            if bit:
                rows[row][column] = ord("1")
                rows[column][row] = ord("1")
            edge_index += 1
    return rows


def matrix_input(record: str) -> bytes:
    rows = decode_graph6(record)
    order = len(rows)
    return (
        f"n={order}\n".encode("ascii")
        + b"\n".join(bytes(row) for row in rows)
        + b"\n"
    )


def read_graph6_inputs(paths: Sequence[str]) -> list[str]:
    requested = list(paths) if paths else ["-"]
    if requested.count("-") > 1:
        raise ExploreError("standard input may be named only once")

    records: list[str] = []
    for path_text in requested:
        if path_text == "-":
            stream = sys.stdin
            display_name = "standard input"
            close_stream = False
        else:
            path = Path(path_text)
            try:
                stream = path.open(encoding="ascii")
            except (OSError, UnicodeError) as error:
                raise ExploreError(f"cannot read '{path_text}': {error}") from error
            display_name = path_text
            close_stream = True

        try:
            for line_number, source_line in enumerate(stream, start=1):
                line = source_line.strip()
                if line.startswith(GRAPH6_HEADER):
                    line = line[len(GRAPH6_HEADER) :]
                if not line:
                    continue
                try:
                    decode_graph6(line)
                except ExploreError as error:
                    raise ExploreError(
                        f"{display_name}:{line_number}: {error}"
                    ) from error
                records.append(line)
        finally:
            if close_stream:
                stream.close()

    if not records:
        raise ExploreError("no graph6 graphs were supplied")
    return records


def canonicalize_graphs(
    records: Sequence[str], labelg: Sequence[str]
) -> list[str]:
    if not records:
        return []
    payload = ("\n".join(records) + "\n").encode("ascii")
    result = subprocess.run(
        [*labelg, "-qg"],
        input=payload,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if result.returncode != 0:
        detail = result.stderr.decode("utf-8", "replace").strip()
        raise ExploreError(
            f"nauty labelg failed with status {result.returncode}"
            + (f": {detail}" if detail else "")
        )
    try:
        canonical = [
            line.decode("ascii")
            for line in result.stdout.splitlines()
            if line and line != GRAPH6_HEADER.encode("ascii")
        ]
    except UnicodeDecodeError as error:
        raise ExploreError("nauty labelg produced non-ASCII output") from error
    if len(canonical) != len(records):
        raise ExploreError(
            "nauty labelg did not return one graph for every input graph"
        )
    for record in canonical:
        decode_graph6(record)
    return sorted(set(canonical))


def prepare_output_directory(path: Path) -> Path:
    try:
        path.mkdir(parents=True, exist_ok=True)
    except OSError as error:
        raise ExploreError(f"cannot create output directory '{path}': {error}")
    try:
        entries = list(path.iterdir())
    except OSError as error:
        raise ExploreError(f"cannot inspect output directory '{path}': {error}")
    if entries:
        raise ExploreError(f"output directory '{path}' is not empty")
    work = path / ".work"
    try:
        work.mkdir()
    except OSError as error:
        raise ExploreError(f"cannot create working directory '{work}': {error}")
    return work


def write_graph6_file(path: Path, records: Iterable[str]) -> None:
    temporary = path.with_name(f".{path.name}.tmp")
    ordered = sorted(records)
    try:
        with temporary.open("w", encoding="ascii", newline="\n") as output:
            for record in ordered:
                output.write(record)
                output.write("\n")
        temporary.replace(path)
    except OSError as error:
        try:
            temporary.unlink()
        except OSError:
            pass
        raise ExploreError(f"cannot write '{path}': {error}") from error


def terminate_processes(processes: Sequence[subprocess.Popen[bytes]]) -> None:
    for process in reversed(processes):
        if process.poll() is None:
            try:
                process.send_signal(signal.SIGTERM)
            except ProcessLookupError:
                pass
    finish_by = time.monotonic() + 1.0
    for process in processes:
        remaining = finish_by - time.monotonic()
        if process.poll() is None and remaining > 0:
            try:
                process.wait(timeout=remaining)
            except subprocess.TimeoutExpired:
                pass
    for process in reversed(processes):
        if process.poll() is None:
            try:
                process.kill()
            except ProcessLookupError:
                pass
    for process in processes:
        try:
            process.wait()
        except ChildProcessError:
            pass


def run_switching_job(
    job_number: int,
    record: str,
    round_number: int,
    label_jobs: int,
    switching: Switching,
    toolchain: Toolchain,
    work_dir: Path,
    stop_event: threading.Event,
) -> JobResult:
    stem = f"round-{round_number:04d}-job-{job_number:08d}"
    output_paths = tuple(
        work_dir / f"{stem}-label-{label_number:04d}.g6"
        for label_number in range(label_jobs)
    )
    error_path = work_dir / f"{stem}.err"
    processes: list[subprocess.Popen[bytes]] = []
    output_files = []
    environment = os.environ.copy()
    environment["LC_ALL"] = "C"
    environment["TMPDIR"] = str(work_dir)

    if stop_event.is_set():
        raise concurrent.futures.CancelledError()

    try:
        with error_path.open("w+b") as errors:
            for output_path in output_paths:
                output_files.append(output_path.open("wb"))
            switch_command = [
                *toolchain.graphswitching,
                *switching.arguments(),
                "--format",
                "graph6",
            ]
            switch = subprocess.Popen(
                switch_command,
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=errors,
                env=environment,
            )
            processes.append(switch)
            assert switch.stdout is not None
            label_inputs = []
            commands = [switch_command]
            sorts = []
            for output in output_files:
                label_command = [*toolchain.labelg, "-qg"]
                label = subprocess.Popen(
                    label_command,
                    stdin=subprocess.PIPE,
                    stdout=subprocess.PIPE,
                    stderr=errors,
                    env=environment,
                )
                processes.append(label)
                commands.append(label_command)
                assert label.stdin is not None
                assert label.stdout is not None
                label_inputs.append(label.stdin)

                sort_command = [*toolchain.sort, "-u"]
                sort = subprocess.Popen(
                    sort_command,
                    stdin=label.stdout,
                    stdout=output,
                    stderr=errors,
                    env=environment,
                )
                processes.append(sort)
                commands.append(sort_command)
                sorts.append(sort)
                label.stdout.close()

            assert switch.stdin is not None
            try:
                switch.stdin.write(matrix_input(record))
                switch.stdin.close()
            except BrokenPipeError:
                pass

            dispatch_errors: list[BaseException] = []

            def dispatch_graphs() -> None:
                try:
                    for graph_number, line in enumerate(switch.stdout):
                        label_inputs[graph_number % label_jobs].write(line)
                except BaseException as error:
                    dispatch_errors.append(error)
                finally:
                    switch.stdout.close()
                    for label_input in label_inputs:
                        try:
                            label_input.close()
                        except OSError:
                            pass

            dispatcher = threading.Thread(target=dispatch_graphs)
            dispatcher.start()
            while dispatcher.is_alive():
                if stop_event.wait(0.1):
                    terminate_processes(processes)
                    dispatcher.join()
                    raise concurrent.futures.CancelledError()
            dispatcher.join()
            if dispatch_errors:
                terminate_processes(processes)
                raise dispatch_errors[0]

            while any(sort.poll() is None for sort in sorts):
                if stop_event.wait(0.1):
                    terminate_processes(processes)
                    raise concurrent.futures.CancelledError()

            statuses = [process.wait() for process in processes]
            if any(status != 0 for status in statuses):
                errors.flush()
                errors.seek(0)
                detail = errors.read().decode("utf-8", "replace").strip()
                failures = [
                    f"{shlex.join(command)} exited {status}"
                    for command, status in zip(commands, statuses)
                    if status != 0
                ]
                raise ExploreError(
                    "; ".join(failures) + (f": {detail}" if detail else "")
                )
    except (OSError, subprocess.SubprocessError) as error:
        terminate_processes(processes)
        raise ExploreError(f"could not run switching pipeline: {error}") from error
    except BaseException:
        terminate_processes(processes)
        for output_path in output_paths:
            try:
                output_path.unlink()
            except OSError:
                pass
        raise
    finally:
        for output in output_files:
            output.close()

    return JobResult(output_paths=output_paths, error_path=error_path)


def dreadnaut_graph(record: str) -> str:
    rows = decode_graph6(record)
    lines = [f"n={len(rows)}", "g"]
    for vertex, row in enumerate(rows):
        neighbours = [
            str(other)
            for other in range(vertex + 1, len(rows))
            if row[other] == ord("1")
        ]
        lines.append(
            f"{vertex}:{' ' if neighbours else ''}{' '.join(neighbours)};"
        )
    lines.append(".")
    lines.append("x")
    return "\n".join(lines)


def automorphism_sizes(
    records: Sequence[str],
    dreadnaut: Sequence[str],
    timeout: float | None,
) -> list[decimal.Decimal]:
    if not records:
        return []
    payload = "-a -m\n" + "\n".join(
        dreadnaut_graph(record) for record in records
    ) + "\nq\n"
    try:
        result = subprocess.run(
            list(dreadnaut),
            input=payload.encode("ascii"),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=timeout,
            check=False,
        )
    except subprocess.TimeoutExpired:
        raise
    except (OSError, subprocess.SubprocessError) as error:
        raise ExploreError(f"could not run nauty dreadnaut: {error}") from error
    if result.returncode != 0:
        detail = (result.stdout + result.stderr).decode(
            "utf-8", "replace"
        ).strip()
        raise ExploreError(
            f"nauty dreadnaut failed with status {result.returncode}"
            + (f": {detail}" if detail else "")
        )
    output = (result.stdout + result.stderr).decode("utf-8", "replace")
    matches = AUTOMORPHISM_SIZE_PATTERN.findall(output)
    if len(matches) != len(records):
        raise ExploreError(
            "could not read one automorphism group size per graph from "
            "nauty dreadnaut"
        )
    try:
        return [decimal.Decimal(value.strip()) for value in matches]
    except decimal.InvalidOperation as error:
        raise ExploreError(
            "nauty dreadnaut returned an invalid automorphism group size"
        ) from error


def remaining_seconds(deadline: float | None) -> float | None:
    if deadline is None:
        return None
    return max(0.0, deadline - time.monotonic())


def apply_automorphism_filter(
    candidates: Sequence[str],
    invariant_filter: AutomorphismFilter | None,
    dreadnaut: Sequence[str] | None,
    deadline: float | None,
) -> tuple[list[str], bool]:
    if invariant_filter is None:
        return list(candidates), False
    if dreadnaut is None:
        raise ExploreError("automorphism filtering requires dreadnaut")

    kept: list[str] = []
    for offset in range(0, len(candidates), AUTOMORPHISM_BATCH_SIZE):
        timeout = remaining_seconds(deadline)
        if timeout is not None and timeout <= 0.0:
            return kept, True
        batch = candidates[offset : offset + AUTOMORPHISM_BATCH_SIZE]
        try:
            sizes = automorphism_sizes(batch, dreadnaut, timeout)
        except subprocess.TimeoutExpired:
            return kept, True
        kept.extend(
            record
            for record, size in zip(batch, sizes)
            if invariant_filter.accepts(size)
        )
    return kept, False


def read_job_records(result: JobResult) -> list[str]:
    records: set[str] = set()
    try:
        for output_path in result.output_paths:
            with output_path.open(encoding="ascii") as source:
                records.update(line.strip() for line in source if line.strip())
    except (OSError, UnicodeError) as error:
        raise ExploreError(
            f"cannot read switching output: {error}"
        ) from error
    finally:
        for output_path in result.output_paths:
            try:
                output_path.unlink()
            except OSError:
                pass
        try:
            result.error_path.unlink()
        except OSError:
            pass
    for record in records:
        decode_graph6(record)
    return sorted(records)


def explore_round(
    round_number: int,
    frontier: Sequence[str],
    seen: set[str],
    jobs: int,
    maximum_graphs: int | None,
    switching: Switching,
    toolchain: Toolchain,
    invariant_filter: AutomorphismFilter | None,
    work_dir: Path,
    deadline: float | None,
) -> RoundResult:
    stop_event = threading.Event()
    new_graphs: set[str] = set()
    processed = 0
    timed_out = False
    graph_limit_reached = False
    next_job = iter(enumerate(frontier))
    active: dict[concurrent.futures.Future[JobResult], int] = {}
    switching_jobs = min(jobs, len(frontier))
    base_label_jobs, extra_label_jobs = divmod(jobs, switching_jobs)

    def submit_one(
        executor: concurrent.futures.ThreadPoolExecutor,
    ) -> bool:
        try:
            job_number, record = next(next_job)
        except StopIteration:
            return False
        future = executor.submit(
            run_switching_job,
            job_number,
            record,
            round_number,
            base_label_jobs + (job_number < extra_label_jobs),
            switching,
            toolchain,
            work_dir,
            stop_event,
        )
        active[future] = job_number
        return True

    with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as executor:
        for _ in range(switching_jobs):
            submit_one(executor)

        try:
            while active:
                timeout = remaining_seconds(deadline)
                if timeout is not None and timeout <= 0.0:
                    timed_out = True
                    stop_event.set()
                    break
                done, _ = concurrent.futures.wait(
                    active,
                    timeout=timeout,
                    return_when=concurrent.futures.FIRST_COMPLETED,
                )
                if not done:
                    timed_out = True
                    stop_event.set()
                    break

                for future in done:
                    active.pop(future)
                    result = future.result()
                    processed += 1
                    candidates = [
                        record
                        for record in read_job_records(result)
                        if record not in seen and record not in new_graphs
                    ]
                    candidates, filter_timed_out = apply_automorphism_filter(
                        candidates,
                        invariant_filter,
                        toolchain.dreadnaut,
                        deadline,
                    )
                    if filter_timed_out:
                        timed_out = True
                        stop_event.set()
                    for record in candidates:
                        if (
                            maximum_graphs is not None
                            and len(seen) + len(new_graphs)
                            >= maximum_graphs
                        ):
                            graph_limit_reached = True
                            stop_event.set()
                            break
                        new_graphs.add(record)
                        if (
                            maximum_graphs is not None
                            and len(seen) + len(new_graphs)
                            >= maximum_graphs
                        ):
                            graph_limit_reached = True
                            stop_event.set()
                            break
                    if timed_out or graph_limit_reached:
                        break
                    submit_one(executor)
                if timed_out or graph_limit_reached:
                    break
        except BaseException:
            stop_event.set()
            raise
        finally:
            if stop_event.is_set():
                for future in active:
                    future.cancel()

    return RoundResult(
        graphs=new_graphs,
        processed=processed,
        timed_out=timed_out,
        graph_limit_reached=graph_limit_reached,
    )


def remove_work_directory(work_dir: Path) -> None:
    try:
        for path in work_dir.iterdir():
            try:
                path.unlink()
            except OSError:
                pass
        work_dir.rmdir()
    except OSError:
        pass


def report(message: str, quiet: bool) -> None:
    if not quiet:
        print(message, file=sys.stderr, flush=True)


def run(options: argparse.Namespace) -> int:
    toolchain = resolve_toolchain(options)
    raw_seeds = read_graph6_inputs(options.inputs)
    seeds = canonicalize_graphs(raw_seeds, toolchain.labelg)
    output_dir = Path(options.output_dir)
    work_dir = prepare_output_directory(output_dir)
    invariant_filter = None
    if any(
        value is not None
        for value in (
            options.aut_size,
            options.min_aut_size,
            options.max_aut_size,
        )
    ):
        invariant_filter = AutomorphismFilter(
            exact=options.aut_size,
            minimum=options.min_aut_size,
            maximum=options.max_aut_size,
        )
    switching = Switching(
        method=options.method,
        part_size=options.part_size,
        symmetry=options.sym,
    )

    try:
        if options.max_graphs is not None and len(seeds) > options.max_graphs:
            raise ExploreError(
                f"{len(seeds)} canonical seeds exceed --max-graphs="
                f"{options.max_graphs}"
            )
        seen = set(seeds)
        frontier = list(seeds)
        write_graph6_file(output_dir / "initial.g6", seeds)
        write_graph6_file(output_dir / "all.g6", seen)
        report(
            f"initial: {len(raw_seeds)} input graph(s), "
            f"{len(seeds)} canonical class(es)",
            options.quiet,
        )

        if options.max_graphs is not None and len(seen) >= options.max_graphs:
            report(
                f"stopped: graph limit reached ({len(seen)} total)",
                options.quiet,
            )
            return 0

        deadline = (
            time.monotonic() + options.time_limit
            if options.time_limit is not None
            else None
        )
        stop_reason = f"completed {options.rounds} round(s)"
        completed_rounds = 0
        for round_number in range(1, options.rounds + 1):
            result = explore_round(
                round_number=round_number,
                frontier=frontier,
                seen=seen,
                jobs=options.jobs,
                maximum_graphs=options.max_graphs,
                switching=switching,
                toolchain=toolchain,
                invariant_filter=invariant_filter,
                work_dir=work_dir,
                deadline=deadline,
            )
            write_graph6_file(
                output_dir / f"round-{round_number:04d}.g6",
                result.graphs,
            )
            seen.update(result.graphs)
            write_graph6_file(output_dir / "all.g6", seen)
            completed_rounds = round_number
            report(
                f"round {round_number}: processed {result.processed}/"
                f"{len(frontier)}, found {len(result.graphs)} new, "
                f"{len(seen)} total",
                options.quiet,
            )
            frontier = sorted(result.graphs)

            if result.timed_out:
                stop_reason = "time limit reached"
                break
            if result.graph_limit_reached:
                stop_reason = "graph limit reached"
                break
            if not frontier:
                stop_reason = "no new graph classes"
                break

        if options.rounds == 0:
            stop_reason = "completed 0 rounds"
        report(
            f"stopped after {completed_rounds} round(s): {stop_reason}; "
            f"{len(seen)} graph class(es) in all.g6",
            options.quiet,
        )
        return 0
    finally:
        remove_work_directory(work_dir)


def main(arguments: Sequence[str] | None = None) -> int:
    try:
        options = parse_arguments(arguments)
        return run(options)
    except ExploreError as error:
        print(f"{program_name()}: {error}", file=sys.stderr)
        return 1
    except KeyboardInterrupt:
        print(f"{program_name()}: interrupted", file=sys.stderr)
        return 130


if __name__ == "__main__":
    raise SystemExit(main())
