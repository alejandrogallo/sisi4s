include $(top_srcdir)/etc/make/emacs.mk

# variables
#
SPHINX_BUILD_DIR = build
DOXYGEN_INDEX = $(SPHINX_BUILD_DIR)/doxygen/html/index.html

DX_EXTRA_STYLESHEET_URL = https://github.com/jothepro/doxygen-awesome-css.git
DX_HTML_EXTRA_FOLDER = $(builddir)/doxygen-awesome-css
DX_HTML_EXTRA_STYLESHEET = $(DX_HTML_EXTRA_FOLDER)/doxygen-awesome.css


all: doxygen sphinx
.PHONY: sphinx doxygen

doxygen: $(DOXYGEN_INDEX)

noinst_PROGRAMS = index.rst manual.rst
index_rst_SOURCES = $(top_srcdir)/README.rst
manual_rst_SOURCES = ./manual.org

$(DOXYGEN_INDEX): Doxyfile $(DX_HTML_EXTRA_STYLESHEET)
	mkdir -p $(@D)
	doxygen $<

$(DX_HTML_EXTRA_STYLESHEET):
	git clone $(DX_EXTRA_STYLESHEET_URL) $(DX_HTML_EXTRA_FOLDER)

$(SPHINX_BUILD_DIR)/index.html: index.rst manual.rst conf.py
	sphinx-build -b html . $(SPHINX_BUILD_DIR)
sphinx: $(SPHINX_BUILD_DIR)/index.html

index.rst: $(top_srcdir)/README.rst
	cp $< $@

$(builddir)/manual.rst: manual.org
	$(EMACS_HTML) $< -f org-rst-export-to-rst
	mv $(dir $<)/manual.rst $@

clean-local:
	rm -rfv sphinx doxygen/ index.rst manual.rst \
		$(SPHINX_BUILD_DIR) \
		$(DX_HTML_EXTRA_FOLDER)
