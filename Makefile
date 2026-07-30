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

VERSION_FILE := VERSION
PROJECT_VERSION := $(shell sed -n '1p' $(VERSION_FILE))
VERSION_CPPFLAGS := -DGRAPHSWITCHING_VERSION=\"$(PROJECT_VERSION)\"
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

graphswitching: $(TOOL_SOURCE) $(LIBRARY_SOURCE) \
		$(INCLUDE_DIR)/graphswitching.h $(VERSION_FILE)
	$(CC) $(CPPFLAGS) $(VERSION_CPPFLAGS) -I$(INCLUDE_DIR) \
		$(CFLAGS) $(WARNFLAGS) \
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
	while read -r name vertices part_size; do \
		case "$$name" in \
			""|\#*) continue ;; \
		esac; \
		printf "Checking %s (n=%s, p=%s)... " \
			"$$name" "$$vertices" "$$part_size"; \
		legacy=$$(./gen_all_srgs_wqh_generic \
			"$$vertices" "$$part_size" \
			< "tests/$$name.matrix" | cksum); \
		current=$$(./graphswitching --method wqh \
			--vertices "$$vertices" --part-size "$$part_size" \
			< "tests/$$name.matrix" | cksum); \
		test "$$legacy" = "$$current"; \
		echo "ok"; \
	done < $(TEST_CASES)
	@set -eu; \
	test "$$(./graphswitching --version)" = \
		"graphswitching $$(sed -n '1p' VERSION)"; \
	./graphswitching --help | grep -q -- '--method=METHOD'; \
	gm_auto=$$(./graphswitching \
		< tests/petersen.matrix | cksum); \
	gm_explicit=$$(./graphswitching --method gm --vertices 10 \
		< tests/petersen.matrix | cksum); \
	test "$$gm_auto" = "$$gm_explicit"; \
	gm_header=$$({ printf 'n=10\n'; cat tests/petersen.matrix; } \
		| ./graphswitching | cksum); \
	test "$$gm_auto" = "$$gm_header"; \
	./graphswitching --input tests/petersen.matrix \
		--output /dev/null; \
	if ./graphswitching 10 2 < tests/petersen.matrix \
		> /dev/null 2>&1; then exit 1; fi; \
	if ./graphswitching --part-size 2 < tests/petersen.matrix \
		> /dev/null 2>&1; then exit 1; fi; \
	printf "Checking default GM against gen_all_srgs (Sp(6,2))... "; \
	gm_legacy=$$(./gen_all_srgs \
		< tests/symplectic-sp6-2.matrix | cksum); \
	gm_current=$$(./graphswitching --vertices 63 \
		< tests/symplectic-sp6-2.matrix | cksum); \
	test "$$gm_legacy" = "$$gm_current"; \
	echo "ok"
	@echo "All builds, CLI checks, and equivalent-output checks passed."

install: $(PROGRAMS)
	$(INSTALL) -d "$(DESTDIR)$(BINDIR)"
	$(INSTALL_PROGRAM) $(PROGRAMS) "$(DESTDIR)$(BINDIR)"

uninstall:
	rm -f $(PROGRAMS:%="$(DESTDIR)$(BINDIR)/%")

clean:
	rm -f $(PROGRAMS)
