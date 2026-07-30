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
PROGRAMS := \
	gen_all_srgs \
	gen_all_srgs_64vts \
	gen_all_srgs_wqh6 \
	gen_all_srgs_wqh_generic

.PHONY: all check clean install uninstall

all: $(PROGRAMS)

$(PROGRAMS): %: $(SOURCE_DIR)/%.c
	$(CC) $(CPPFLAGS) $(CFLAGS) $(WARNFLAGS) $< $(LDFLAGS) $(LDLIBS) -o $@

check: $(PROGRAMS)
	@for program in $(PROGRAMS); do \
		test -x "$$program" || exit 1; \
	done
	@echo "All programs built successfully."

install: $(PROGRAMS)
	$(INSTALL) -d "$(DESTDIR)$(BINDIR)"
	$(INSTALL_PROGRAM) $(PROGRAMS) "$(DESTDIR)$(BINDIR)"

uninstall:
	rm -f $(PROGRAMS:%="$(DESTDIR)$(BINDIR)/%")

clean:
	rm -f $(PROGRAMS)
