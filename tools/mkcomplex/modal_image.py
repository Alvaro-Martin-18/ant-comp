"""Modal image for mkcomplex (used by ``--executor modal``).

The task script is pure shell string work -- bash plus one awk one-liner for the
FASTA -- so the image is a bare slim Debian: no scientific stack, no GPU, nothing
pip-installed, and it builds in seconds. The manifest builder sets
``gpus_per_task = 0``.
"""

import modal

RESOURCES = {"cpu": 1, "memory": "1G", "timeout": "00:10:00"}


def image() -> modal.Image:
    return modal.Image.debian_slim(python_version="3.12")
