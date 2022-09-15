EMACS = emacs -q --batch --load $(top_srcdir)/config.el
EMACS_HTML = $(EMACS) --load $(top_srcdir)/etc/emacs/html.el

define tangle
$(EMACS) $(1) --eval '(org-babel-tangle)'
endef

SUFFIXES = .html .org .rst

.org.html:
	$(EMACS_HTML) $< -f org-html-export-to-html

.org.rst:
	$(EMACS_HTML) $< -f org-rst-export-to-rst
