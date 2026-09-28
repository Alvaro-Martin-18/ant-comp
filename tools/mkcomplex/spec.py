from pathlib import Path

from prosapia.core import Tool

from .collect_mkcomplex import collect_mkcomplex
from .run_mkcomplex import add_run_mkcomplex_args, build_mkcomplex_manifest

TOOL = Tool(
    name="mkcomplex",
    action="update",
    description="Rebuild a design's multi-chain complex sequence around fixed chains.",
    default_script=str(Path(__file__).parent / "mkcomplex.sh"),
    default_input_column="proteinmpnn_sequence",
    build_manifest_fn=build_mkcomplex_manifest,
    add_run_args_fn=add_run_mkcomplex_args,
    collect_fn=collect_mkcomplex,
)
