#!/bin/bash
~/miniconda3/envs/gbp1/bin/pip install torch --index-url https://download.pytorch.org/whl/cpu -q 2>&1 | tail -3
~/miniconda3/envs/gbp1/bin/pip install scdiffeq -q 2>&1 | tail -5
~/miniconda3/envs/gbp1/bin/python -c "import scdiffeq; print(\"scdiffeq OK\", scdiffeq.__version__)" 2>&1
