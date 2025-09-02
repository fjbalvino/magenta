SHELL := /bin/bash
PYTHON := python

# Directories
DATA_DIR ?= data
RESULTS_DIR ?= results
LOGS_DIR ?= logs

.PHONY: help setup fetch_mangrove fetch_non_mangrove download_mangrove download_non_mangrove fastqc assemble qc all

help:
	@echo "Targets:"
	@echo "  setup                Create folders"
	@echo "  fetch_mangrove       Buscar metadatos de manglar"
	@echo "  fetch_non_mangrove   Buscar metadatos de no-manglar"
	@echo "  download_mangrove    Descargar/convertir datos de manglar"
	@echo "  download_non_mangroveDescargar/convertir datos de no-manglar"
	@echo "  fastqc               Correr FastQC en paralelo"
	@echo "  assemble             Ensamblar (MEGAHIT o MetaSPAdes)"
	@echo "  qc                   (placeholder) MultiQC si está disponible"
	@echo "  all                  Flujo completo"

setup:
	mkdir -p $(DATA_DIR) $(RESULTS_DIR) $(LOGS_DIR)

fetch_mangrove: setup
	$(PYTHON) scripts/01_magenta_fetch_mangrove.py

fetch_non_mangrove: setup
	$(PYTHON) scripts/magenta_fetch_non_mangrove.py

download_mangrove: setup
	$(PYTHON) scripts/02_descargar_y_convertir_mangrove.py

download_non_mangrove: setup
	$(PYTHON) scripts/descargar_y_convertir_no_mangrove.py

fastqc: setup
	$(PYTHON) scripts/04_fastqc_parallel.py

assemble: setup
	$(PYTHON) scripts/06_assembly_serial.py

qc:
	@echo "Agrega tu comando de MultiQC aquí si lo deseas."

all: fetch_mangrove download_mangrove fastqc assemble qc
