SUBDIRS = extern src integration-tests

if WITH_DOCS
SUBDIRS += docs
endif


.PHONY: sisi4s extern

extern:
	$(MAKE) -C extern

sisi4s: src/sisi4s
src/sisi4s:
	$(MAKE) -C src

