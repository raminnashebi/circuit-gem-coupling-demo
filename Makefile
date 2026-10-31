PYTHON ?= python

.PHONY: all results docs test

all: results docs

results:
	$(PYTHON) src/coupling_demo.py

docs:
	$(PYTHON) src/build_report.py
	$(PYTHON) src/build_math.py

test:
	$(PYTHON) -m pytest -q
