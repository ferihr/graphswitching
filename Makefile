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

.PHONY: all check clean install uninstall

all: $(PROGRAMS)

$(LEGACY_PROGRAMS): %: $(SOURCE_DIR)/%.c
	$(CC) $(CPPFLAGS) $(CFLAGS) $(WARNFLAGS) $< $(LDFLAGS) $(LDLIBS) -o $@

graphswitching: $(TOOL_SOURCE) $(LIBRARY_SOURCE) $(INCLUDE_DIR)/graphswitching.h
	$(CC) $(CPPFLAGS) -I$(INCLUDE_DIR) $(CFLAGS) $(WARNFLAGS) \
		$(TOOL_SOURCE) $(LIBRARY_SOURCE) $(LDFLAGS) $(LDLIBS) -o $@

check: $(PROGRAMS)
	@for program in $(PROGRAMS); do \
		test -x "$$program" || exit 1; \
	done
	@set -eu; \
	tmpdir=$$(mktemp -d ./.graphswitching-check.XXXXXX); \
	trap 'rm -rf "$$tmpdir"' EXIT HUP INT TERM; \
	./gen_all_srgs_wqh_generic 4 1 < tests/cycle4.matrix \
		> "$$tmpdir/legacy.out"; \
	./graphswitching 4 1 < tests/cycle4.matrix \
		> "$$tmpdir/graphswitching.out"; \
	cmp "$$tmpdir/legacy.out" "$$tmpdir/graphswitching.out"
	@echo "All programs built and equivalent output verified."

install: $(PROGRAMS)
	$(INSTALL) -d "$(DESTDIR)$(BINDIR)"
	$(INSTALL_PROGRAM) $(PROGRAMS) "$(DESTDIR)$(BINDIR)"

uninstall:
	rm -f $(PROGRAMS:%="$(DESTDIR)$(BINDIR)/%")

clean:
	rm -f $(PROGRAMS)
