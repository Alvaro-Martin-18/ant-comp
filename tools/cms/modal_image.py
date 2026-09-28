"""Modal image for cms (used by ``--executor modal``).

cms-cuda JIT-compiles its kernels with NVRTC through CuPy, which also needs the CUDA
toolkit HEADERS -- the *runtime* image fails every design with "Failed to find CUDA
headers" (verified). So the image starts from NVIDIA's CUDA 12 *devel* image and uses
``cupy-cuda12x``: CUDA 12 runs on every driver Modal currently offers, whereas the ``[gpu]`` extra's
``cupy-cuda13x`` needs a CUDA 13 driver.

cms-cuda is NOT pip-installed: its source dir is ``cms-cuda/`` (hyphen) while its
pyproject looks for ``cms_cuda*``, so ``pip install git+...`` installs metadata and
no importable module. The pinned checkout's package dir is copied onto the path as
``cms_cuda`` instead, and its runtime deps are installed by hand. Revisit if upstream
fixes the packaging.

No weights, no Volume. The first GPU call in a container compiles the kernels
(~60 s on an L4, measured; cached in ~/.cupy/kernel_cache for the container's
lifetime), which is why
the builder packs many designs per task.
"""

import modal

CMS_CUDA_COMMIT = "bbda09cb05a7072dc76b81ad2a01288fe5e63f2b"  # 2026-09-26, v0.2.0
RESOURCES = {"gpu": "L4", "cpu": 4, "memory": "16G", "timeout": "01:00:00"}


def image() -> modal.Image:
    return (
        modal.Image.from_registry(
            "nvidia/cuda:12.6.3-devel-ubuntu22.04", add_python="3.12"
        )
        .apt_install("git")
        .pip_install("numpy", "scipy", "biopython", "cupy-cuda12x")
        .run_commands(
            "git clone https://github.com/ullahsamee/cms-cuda /opt/cms-cuda-src",
            f"git -C /opt/cms-cuda-src checkout {CMS_CUDA_COMMIT}",
            "mkdir -p /opt/cms && cp -r /opt/cms-cuda-src/cms-cuda /opt/cms/cms_cuda",
            "PYTHONPATH=/opt/cms python -c 'import cms_cuda; print(cms_cuda.__version__)'",
        )
        .env({"PYTHONPATH": "/opt/cms"})
    )
