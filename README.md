# 📖 Documentation: https://euclid-surface-photometry.readthedocs.io/en/latest/index.html
# Euclid Surface Photometry Pipeline    

This repository contains the code and configuration files for running the Euclid surface photometry pipeline, as described in *Euclid Consortium: Sanchez-Alarcon et al. (2026)*.

## Overview

The pipeline is designed to process astronomical images and catalogs from the Euclid mission, performing tasks such as cutout creation, profile measurement, and photometric analysis. It is modular, allowing users to enable or disable specific processing steps.

## Features

- Modular pipeline with configurable steps
- Support for custom catalogs and data paths
- Flexible parameterization via YAML configuration
- Logging and output management

## Getting Started

### Prerequisites

- Python 3.8+
- Required Python packages (see `requirements.txt` if available)

### Configuration

Edit the `euclid_config.yaml` file to set up your pipeline steps, input catalogs, data paths, and parameters. Comments in the YAML file explain each option.

You can also generate or modify the configuration file programmatically:

```python
from write_euclid_config import write_euclid_config_yaml

# Write default config
write_euclid_config_yaml("euclid_config.yaml")

# Or override specific values
custom_config = {
    "parameters": {"n_kron": 5, "max_size": 1500},
    "logging": {"level": "DEBUG"}
}
write_euclid_config_yaml("euclid_config.yaml", config_dict=custom_config)