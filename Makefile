# Makefile

# Variables
URL = https://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/mitochondrion.1.genomic.gbff.gz
DATA_DIR = data/dir
VENV_DIR = .venv
REQUIREMENTS_FILE = requirements.txt

# Default target
all:  create_venv run-genome-process

# Create Python virtual environment
create_venv: $(VENV_DIR)/bin/activate

$(VENV_DIR)/bin/activate: $(REQUIREMENTS_FILE)
	python -m venv $(VENV_DIR)
	@echo "Virtual environment created with Python $(PY_VERSION)"
	@echo "Activating virtual environment and installing dependencies..."
	./$(VENV_DIR)/bin/pip install -U pip
	./$(VENV_DIR)/bin/pip install -r $(REQUIREMENTS_FILE)
	@echo "Virtual environment setup complete."

# Run Python script
run-genome-process: venv
	@echo "Running gb_read_utils.py..."
	./$(VENV_DIR)/bin/python src/gb_read_utils.py

# Clean up data directory
clean_data:
	@rm -rf $(DATA_DIR)/*

# Clean up virtual environment
clean_venv:
	@rm -rf $(VENV_DIR)

# Full cleanup
clean: clean_data clean_venv

.PHONY: all  create_venv run-genome-process clean_data clean_venv clean
