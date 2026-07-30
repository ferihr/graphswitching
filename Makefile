CC ?= cc
CPPFLAGS ?=
CFLAGS ?= -O3 -funroll-loops
WARNFLAGS ?= -Wall -Wextra
LDFLAGS ?=
LDLIBS ?=

PREFIX ?= /usr/local
BINDIR ?= $(PREFIX)/bin
INSTALL ?= install
INSTALL_PROGRAM ?= $(INSTALL) -m 0755
PYTHON ?= python3

BENCHMARK_SUITE ?= quick
BENCHMARK_RUNS ?= 1
BENCHMARK_CASES ?=

SOURCE_DIR := code
INCLUDE_DIR := include
LIBRARY_SOURCE := src/graphswitching.c
TOOL_SOURCE := tools/graphswitching.c
LEGACY_PROGRAMS := \
	gen_all_srgs \
	gen_all_srgs_64vts \
	gen_all_srgs_wqh6 \
	gen_all_srgs_wqh_generic
PROGRAMS := graphswitching $(LEGACY_PROGRAMS)
TEST_CASES := tests/cases.txt
BENCHMARK_PROGRAMS := \
	graphswitching \
	gen_all_srgs \
	gen_all_srgs_64vts \
	gen_all_srgs_wqh_generic

.PHONY: all benchmark check clean install uninstall

all: $(PROGRAMS)

$(LEGACY_PROGRAMS): %: $(SOURCE_DIR)/%.c
	$(CC) $(CPPFLAGS) $(CFLAGS) $(WARNFLAGS) $< $(LDFLAGS) $(LDLIBS) -o $@

graphswitching: $(TOOL_SOURCE) $(LIBRARY_SOURCE) $(INCLUDE_DIR)/graphswitching.h
	$(CC) $(CPPFLAGS) -I$(INCLUDE_DIR) $(CFLAGS) $(WARNFLAGS) \
		$(TOOL_SOURCE) $(LIBRARY_SOURCE) $(LDFLAGS) $(LDLIBS) -o $@

benchmark: $(BENCHMARK_PROGRAMS)
	$(PYTHON) benchmarks/run.py \
		--suite "$(BENCHMARK_SUITE)" \
		--runs "$(BENCHMARK_RUNS)" $(BENCHMARK_CASES)

check: $(PROGRAMS) $(TEST_CASES)
	@for program in $(PROGRAMS); do \
		test -x "$$program" || exit 1; \
	done
	@set -eu; \
	tmpdir=$$(mktemp -d ./.graphswitching-check.XXXXXX); \
	trap 'rm -rf "$$tmpdir"' EXIT HUP INT TERM; \
	while read -r name vertices part_size; do \
		case "$$name" in \
			""|\#*) continue ;; \
		esac; \
		printf "Checking %s (n=%s, p=%s)... " \
			"$$name" "$$vertices" "$$part_size"; \
		./gen_all_srgs_wqh_generic "$$vertices" "$$part_size" \
			< "tests/$$name.matrix" > "$$tmpdir/$$name.legacy"; \
		./graphswitching "$$vertices" "$$part_size" \
			< "tests/$$name.matrix" > "$$tmpdir/$$name.graphswitching"; \
		cmp "$$tmpdir/$$name.legacy" \
			"$$tmpdir/$$name.graphswitching"; \
		echo "ok"; \
	done < $(TEST_CASES)
	@echo "All programs built and all equivalent-output checks passed."

install: $(PROGRAMS)
	$(INSTALL) -d "$(DESTDIR)$(BINDIR)"
	$(INSTALL_PROGRAM) $(PROGRAMS) "$(DESTDIR)$(BINDIR)"

uninstall:
	rm -f $(PROGRAMS:%="$(DESTDIR)$(BINDIR)/%")

clean:
	rm -f $(PROGRAMS)
