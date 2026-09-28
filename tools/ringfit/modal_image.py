"""Modal image for ringfit (used by ``--executor modal``).

The worker is pure Python -- gemmi for structure parsing, numpy for the Kabsch
fit, the Shrake-Rupley SASA and the neighbour searches -- so the image is a slim
Debian with those two wheels and nothing else: no GPU, no torch, and it builds
in seconds. The manifest builder sets ``gpus_per_task = 0``.
"""

import modal

RESOURCES = {"cpu": 2, "memory": "8G", "timeout": "00:30:00"}


def image() -> modal.Image:
    return modal.Image.debian_slim(python_version="3.12").pip_install("gemmi", "numpy")
