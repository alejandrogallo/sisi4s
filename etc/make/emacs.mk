EMACS = emacs -q --batch --load config.el
EMACS_HTML = $(EMACS)

define tangle
$(EMACS) $(1) --eval '(org-babel-tangle)'
endef

SUFFIXES = .html .org .rst

.org.html:
	$(EMACS_HTML) $< -f org-html-export-to-html

.org.rst:
	$(EMACS_HTML) $< -f org-rst-export-to-rst
