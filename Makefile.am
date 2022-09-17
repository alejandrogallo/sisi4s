include SiSi4S.mk

SUBDIRS = src

if WITH_DOCS
SUBDIRS += docs
endif

SUBDIRS += integration-tests

.PHONY: sisi4s

sisi4s: src/sisi4s
src/sisi4s:
	$(MAKE) -C src

