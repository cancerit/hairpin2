#!/bin/bash

if [ -z ${TEST_DIR} ]; then
    echo "TEST_DIR not set!"
    exit 1
fi
PKG_DIR=$(python -c "import os;import hairpin2;import inspect;print(os.path.dirname(inspect.getfile(hairpin2)))")

echo "$(python --version)"
echo "Package source directory: ${PKG_DIR}"

pip install \
    pytest==8.2.2 \
    pytest-cov==5.0.0 \
    faker-biology==0.6.4 \
    factory-boy==3.3.1 \
    pysam==0.22 && \
pytest -m "validate" --cov="${PKG_DIR}" "${TEST_DIR}"

