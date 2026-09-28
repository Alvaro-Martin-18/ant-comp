"""Modal image for chainsel (used by ``--executor modal``).

The worker is pure Python -- gemmi to read PDB/mmCIF, select chains and write the
subset back out -- so the image is a slim Debian with that one wheel and nothing
else: no numpy, no GPU, and it builds in seconds. The manifest builder sets
``gpus_per_task = 0``.
"""

import modal

RESOURCES = {"cpu": 1, "memory": "4G", "timeout": "00:20:00"}


def image() -> modal.Image:
    return modal.Image.debian_slim(python_version="3.12").pip_install("gemmi")
