CC ?= cc
CPPFLAGS ?=
CFLAGS ?= -O3 -funroll-loops
WARNFLAGS ?= -Wall -Wextra
LDFLAGS ?=
LDLIBS ?=
PKG_CONFIG ?= pkg-config
WITH_NAUTY ?= auto
NAUTY_PKG_CPPFLAGS := $(shell $(PKG_CONFIG) --cflags nauty 2>/dev/null)
NAUTY_PKG_LDLIBS := $(shell $(PKG_CONFIG) --libs nauty 2>/dev/null)
NAUTY_DEFAULT_HEADER := $(shell \
	printf '%s\n' '#include <nauty.h>' '#include <naugroup.h>' \
		'#include <schreier.h>' \
	| $(CC) $(CPPFLAGS) -x c -E - >/dev/null 2>&1 \
	&& printf 1 || printf 0)
NAUTY_HEADER_DIRS ?= /usr/local/include/nauty /usr/include/nauty
NAUTY_SYSTEM_HEADER_DIR := $(firstword $(foreach dir,$(NAUTY_HEADER_DIRS),\
	$(if $(and $(wildcard $(dir)/nauty.h),\
		$(wildcard $(dir)/naugroup.h),\
		$(wildcard $(dir)/schreier.h)),$(dir))))
NAUTY_SYSTEM_CPPFLAGS := $(if $(NAUTY_SYSTEM_HEADER_DIR),\
	-isystem $(NAUTY_SYSTEM_HEADER_DIR))
NAUTY_CPPFLAGS ?= $(if $(NAUTY_PKG_CPPFLAGS),\
	$(NAUTY_PKG_CPPFLAGS),$(if $(filter 1,$(NAUTY_DEFAULT_HEADER)),,\
	$(NAUTY_SYSTEM_CPPFLAGS)))
NAUTY_LDLIBS ?= $(if $(NAUTY_PKG_LDLIBS),$(NAUTY_PKG_LDLIBS),-lnauty)

ifeq ($(WITH_NAUTY),1)
NAUTY_ENABLED := 1
else ifeq ($(WITH_NAUTY),auto)
NAUTY_ENABLED := $(shell \
	printf '%s\n' '#include <nauty.h>' '#include <naugroup.h>' \
		'#include <schreier.h>' \
		'int main(void) {' \
		'  schreier *group = 0; permnode *generators = 0;' \
		'  grouprec *record = groupptr(FALSE); (void)record;' \
		'  nauty_check(WORDSIZE, 1, 1, NAUTYVERSIONID);' \
		'  newgroup(&group, &generators, 1);' \
		'  freeschreier(&group, &generators); return 0;' \
		'}' \
	| $(CC) $(CPPFLAGS) $(NAUTY_CPPFLAGS) -x c - -x none \
		$(LDFLAGS) $(NAUTY_LDLIBS) -o /dev/null \
		>/dev/null 2>&1 && printf 1 || printf 0)
else
NAUTY_ENABLED := 0
endif

ifeq ($(NAUTY_ENABLED),1)
GRAPHSWITCHING_CPPFLAGS := \
	$(NAUTY_CPPFLAGS) -DGRAPHSWITCHING_WITH_NAUTY
GRAPHSWITCHING_LDLIBS := $(NAUTY_LDLIBS)
endif

PREFIX ?= /usr/local
BINDIR ?= $(PREFIX)/bin
INSTALL ?= install
INSTALL_PROGRAM ?= $(INSTALL) -m 0755
PYTHON ?= python3

BENCHMARK_SUITE ?= short
BENCHMARK_RUNS ?= 1
BENCHMARK_CASES ?=
TUNE_GCC_ARGS ?=

VERSION_FILE := VERSION
PROJECT_VERSION := $(shell sed -n '1p' $(VERSION_FILE))
VERSION_CPPFLAGS := -DGRAPHSWITCHING_VERSION=\"$(PROJECT_VERSION)\"
SOURCE_DIR := code
INCLUDE_DIR := include
LIBRARY_SOURCES := \
	src/graphswitching.c \
	src/fixed.c \
	src/gm.c \
	src/wqh.c \
	src/symmetry.c \
	src/io.c \
	src/switching_methods.c
INTERNAL_HEADERS := \
	src/graphswitching_internal.h \
	src/switching_methods.h
TOOL_SOURCE := tools/graphswitching.c
EXPLORE_SOURCE := tools/graphswitching_explore.py
LEGACY_PROGRAMS := \
	gen_all_srgs \
	gen_all_srgs_64vts \
	gen_all_srgs_wqh6 \
	gen_all_srgs_wqh_generic
PROGRAMS := graphswitching graphswitching-explore $(LEGACY_PROGRAMS)
TEST_CASES := tests/cases.txt
BENCHMARK_PROGRAMS := \
	graphswitching \
	gen_all_srgs \
	gen_all_srgs_64vts \
	gen_all_srgs_wqh_generic

.PHONY: all benchmark check check-symmetry clean install tune-gcc uninstall FORCE

all: $(PROGRAMS)

$(LEGACY_PROGRAMS): %: $(SOURCE_DIR)/%.c
	$(CC) $(CPPFLAGS) $(CFLAGS) $(WARNFLAGS) $< $(LDFLAGS) $(LDLIBS) -o $@

graphswitching: FORCE $(TOOL_SOURCE) $(LIBRARY_SOURCES) \
		$(INCLUDE_DIR)/graphswitching.h $(INTERNAL_HEADERS) $(VERSION_FILE)
	$(CC) $(CPPFLAGS) $(GRAPHSWITCHING_CPPFLAGS) \
		$(VERSION_CPPFLAGS) -I$(INCLUDE_DIR) \
		$(CFLAGS) $(WARNFLAGS) \
		$(TOOL_SOURCE) $(LIBRARY_SOURCES) $(LDFLAGS) \
		$(LDLIBS) $(GRAPHSWITCHING_LDLIBS) -o $@

graphswitching-explore: $(EXPLORE_SOURCE) $(VERSION_FILE)
	sed 's/@GRAPHSWITCHING_VERSION@/$(PROJECT_VERSION)/g' $< > $@
	chmod +x $@

FORCE:

benchmark: $(BENCHMARK_PROGRAMS)
	$(PYTHON) benchmarks/run.py \
		--suite "$(BENCHMARK_SUITE)" \
		--runs "$(BENCHMARK_RUNS)" $(BENCHMARK_CASES)

tune-gcc:
	$(PYTHON) benchmarks/tune_gcc.py $(TUNE_GCC_ARGS)

check-symmetry: graphswitching
	@if test "$(NAUTY_ENABLED)" != 1; then \
		echo "check-symmetry: graphswitching was built without nauty" \
			>&2; \
		exit 2; \
	fi
	$(PYTHON) benchmarks/run.py \
		sp6-gm-sym \
		petersen-gm-sym \
		clebsch-gm-sym \
		gq2-4-gm-sym \
		sp6-wqh-sym-p2 \
		gq2-4-wqh-sym-p5 \
		bil223-wqh-sym-p3 \
		sp4-4-wqh-sym-p4 \
		sp6-gm2-aut8-gm-sym \
		sp6-gm3-aut1-gm-sym
	$(PYTHON) tests/check_fixed_methods.py

check: $(PROGRAMS) $(TEST_CASES)
	$(PYTHON) tests/generate_algebraic_fixtures.py
	$(PYTHON) tests/check_fixed_methods.py --data-only
	$(PYTHON) tests/check_graph6_output.py
	PYTHONDONTWRITEBYTECODE=1 \
		$(PYTHON) tests/test_graphswitching_explore.py
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
	help=$$(./graphswitching --help); \
	for option in --method --vertices --part-size --format --input --output \
		--sym --help --version; do \
		printf '%s\n' "$$help" | grep -q -- "$$option"; \
	done; \
	explore_help=$$(./graphswitching-explore --help); \
	for option in --output-dir --jobs --processes --rounds --time-limit \
		--max-graphs --method --part-size --sym --aut-size \
		--min-aut-size --max-aut-size --help --version; do \
		printf '%s\n' "$$explore_help" | grep -q -- "$$option"; \
	done; \
	test "$$(./graphswitching-explore --version)" = \
		"graphswitching-explore $$(sed -n '1p' VERSION)"; \
	printf '%s\n' "$$help" | grep -q -- 'P,P,N-2P'; \
	printf '%s\n' "$$help" | grep -q -- '2P <= N'; \
	printf '%s\n' "$$help" | grep -q -- 'Exit status:'; \
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
	if ./graphswitching --method fano --part-size 2 \
		< tests/petersen.matrix \
		> /dev/null 2>&1; then exit 1; fi; \
	for method in gm wqh ah gm2 is3 is5 fano; do \
		printf '%s\n' "$$help" | grep -q -- "$$method"; \
	done; \
	if test "$(NAUTY_ENABLED)" = 1; then \
		./graphswitching --method gm --sym \
			< tests/petersen.matrix > /dev/null; \
		./graphswitching --method wqh --sym \
			< tests/petersen.matrix > /dev/null; \
		auto=$$(./graphswitching --method gm --vertices 63 \
			--sym=auto --format graph6 \
			< tests/symplectic-sp6-2.matrix | cksum); \
		forced=$$(./graphswitching --method gm --vertices 63 \
			--sym --format graph6 \
			< tests/symplectic-sp6-2.matrix | cksum); \
		test "$$auto" = "$$forced"; \
		auto=$$(./graphswitching --method gm --vertices 63 \
			--sym=auto --format graph6 \
			< tests/symplectic-sp6-2-gm2-aut8.matrix | cksum); \
		forced=$$(./graphswitching --method gm --vertices 63 \
			--sym --format graph6 \
			< tests/symplectic-sp6-2-gm2-aut8.matrix | cksum); \
		test "$$auto" = "$$forced"; \
		auto=$$(./graphswitching --method gm --vertices 63 \
			--sym=auto --format graph6 \
			< tests/symplectic-sp6-2-gm3-aut1.matrix | cksum); \
		ordinary=$$(./graphswitching --method gm --vertices 63 \
			--format graph6 \
			< tests/symplectic-sp6-2-gm3-aut1.matrix | cksum); \
		test "$$auto" = "$$ordinary"; \
	else \
		if ./graphswitching --method gm --sym \
			< tests/petersen.matrix > /dev/null 2>&1; then \
			exit 1; \
		fi; \
		if ./graphswitching --method wqh --sym \
			< tests/petersen.matrix > /dev/null 2>&1; then \
			exit 1; \
		fi; \
		if ./graphswitching --method gm --sym=auto \
			< tests/petersen.matrix > /dev/null 2>&1; then \
			exit 1; \
		fi; \
	fi; \
	if ./graphswitching --vertices 968 \
		< tests/petersen.matrix > /dev/null 2>&1; then exit 1; fi; \
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
