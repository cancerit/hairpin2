# hairpin-core

Maintainable, transparent, implementation of the hairpin detection and flagging algorithm concieved by Mathijs' Sanders. Implemented by Peter Campbell and Alex Byrne

### REQUIREMENTS

* Python3 - tested with 3.12

### INSTALLATION


Within a virtual environment:
```
python -m venv .env
source .env/bin/activate
pip install .
hairpin -h
```

For system-wide access:
```
export INST_PATH=/path/to/install/location/
mkdir -p $INST_PATH
pip install . --target $INST_PATH
export PATH=${PATH}:${INST_PATH}/bin
hairpin -h
```
