# Benchmarks

These benchmarks measure the complete workflow used for switching
experiments:

```text
switching program | amtog -q | labelg -gqt | sort -u | wc -l
```

The runner checks the final number of pairwise nonisomorphic graphs and
reports elapsed wall time plus the total user and system CPU time used by
the pipeline. With multiple repetitions, it reports the median of each
measurement.

## Requirements

- Python 3.8 or newer
- nauty and Traces command-line tools
- the standard `sort` and `wc` utilities

Build and run the quick suite:

```sh
make benchmark
```

The runner automatically looks for both common nauty naming schemes:

```text
amtog, labelg
nauty-amtog, nauty-labelg
```

To choose a prefix explicitly:

```sh
make benchmark NAUTY_PREFIX=nauty-
python3 benchmarks/run.py --nauty-prefix nauty-
```

For paths or nonstandard names, set both command overrides:

```sh
AMTOG=/opt/nauty/bin/amtog \
LABELG=/opt/nauty/bin/labelg \
make benchmark
```

## Suites and cases

The default `quick` suite compares the two original GM programs, the default
GM method in `graphswitching`, and both generic WQH implementations on
Sp(6,2). This reproduces and extends the end-to-end benchmark in the project
discussion.

The `full` suite also uses `p=3` on Sp(6,2) and Paley(73). Those cases emit
no switched matrices, so their running time is dominated by candidate-tuple
enumeration instead of matrix conversion or canonical labeling.

```sh
make benchmark BENCHMARK_SUITE=full
make benchmark BENCHMARK_SUITE=full BENCHMARK_RUNS=3
```

List the cases or run only named cases:

```sh
python3 benchmarks/run.py --list
make benchmark BENCHMARK_CASES="sp6-wqh-p2 paley73-wqh-p3"
python3 benchmarks/run.py --runs 5 sp6-gm sp6-gm-packed
```

The case definitions live in `cases.tsv`. The `quick` suite currently
expects two isomorphism classes from each Sp(6,2) benchmark. The full
suite's `p=3` cases expect zero, matching their intentional
tuple-enumeration-only role.

For useful comparisons, run benchmarks on an otherwise idle machine and
record the compiler, `CFLAGS`, nauty version, CPU, and operating system.
Wall-clock timings from different machines or builds should not be compared
without that context.
