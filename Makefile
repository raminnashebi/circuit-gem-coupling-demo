PYTHON ?= python

.PHONY: all results docs slides method test

all: results docs slides method

results:
	$(PYTHON) src/coupling_demo.py

docs:
	$(PYTHON) src/build_report.py
	$(PYTHON) src/build_math.py

slides:
	$(PYTHON) src/slide_figures.py

method:
	$(PYTHON) src/method_figures.py

test:
	$(PYTHON) -m pytest -q
