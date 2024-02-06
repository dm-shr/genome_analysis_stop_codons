# Makefile

# Variables
URL = https://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/mitochondrion.1.genomic.gbff.gz
DATA_DIR = data/dir
VENV_DIR = .venv
REQUIREMENTS_FILE = requirements.txt

# Default target
all:  create_venv run-genome-process
# download unpack
# # Download the file
# download:
# 	@mkdir -p $(DATA_DIR)
# 	@wget -c $(URL) -O $(DATA_DIR)/mitochondrion.1.genomic.gbff.gz

# # Unpack the .gz file and remove the compressed file
# unpack:
# 	@gunzip -f $(DATA_DIR)/mitochondrion.1.genomic.gbff.gz

# Create Python virtual environment
create_venv: $(VENV_DIR)/bin/activate

$(VENV_DIR)/bin/activate: $(REQUIREMENTS_FILE)
	python -m venv $(VENV_DIR)
	@echo "Virtual environment created with Python $(PY_VERSION)"
	@echo "Activating virtual environment and installing dependencies..."
	./$(VENV_DIR)/bin/pip install -U pip
	./$(VENV_DIR)/bin/pip install -r $(REQUIREMENTS_FILE)
	@echo "Virtual environment setup complete."

# # Install Python requirements
# install_requirements:
# 	@$(VENV_DIR)/bin/pip install -r $(REQUIREMENTS_FILE)

# Run Python script
run-genome-process: venv
	@echo "Running genome_process.py..."
	./$(VENV_DIR)/bin/python genome_process.py

# Clean up data directory
clean_data:
	@rm -rf $(DATA_DIR)/*

# Clean up virtual environment
clean_venv:
	@rm -rf $(VENV_DIR)

# Full cleanup
clean: clean_data clean_venv
# download unpack
.PHONY: all  create_venv run-genome-process clean_data clean_venv clean
