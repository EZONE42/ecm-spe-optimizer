.PHONY: help install install-dev test run optimize transference report all clean

PYTHON ?= python

help:
	@echo "ECM-SPE: Solid Polymer Electrolyte Composition Optimization"
	@echo ""
	@echo "Targets:"
	@echo "  install       Install package in editable mode"
	@echo "  install-dev   Install package with development/test dependencies"
	@echo "  test          Run unit and regression test suite with pytest"
	@echo "  optimize      Run 4-resistance ECM composition optimization sweep"
	@echo "  transference  Run Bruce-Vincent transference number calculation"
	@echo "  report        Print cell 4-resistance breakdown & bottleneck analysis"
	@echo "  all           Run all pipeline steps"
	@echo "  clean         Remove build artifacts and cache files"

install:
	$(PYTHON) -m pip install -e .

install-dev:
	$(PYTHON) -m pip install -e ".[dev]"

test:
	pytest -v

optimize:
	$(PYTHON) -m ecm_spe optimize

transference:
	$(PYTHON) -m ecm_spe transference

report:
	$(PYTHON) -m ecm_spe report

all:
	$(PYTHON) -m ecm_spe all

clean:
	find . -type d -name "__pycache__" -exec rm -rf {} +
	find . -type f -name "*.pyc" -delete
	rm -rf build/ dist/ *.egg-info .pytest_cache/
