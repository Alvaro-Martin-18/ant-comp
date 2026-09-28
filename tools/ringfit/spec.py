from pathlib import Path

from prosapia.core import Tool

from .collect_ringfit import collect_ringfit
from .run_ringfit import add_run_ringfit_args, build_ringfit_manifest

TOOL = Tool(
    name="ringfit",
    action="update",
    description="Score how a two-protomer binder fits the full oligomeric assembly.",
    default_script=str(Path(__file__).parent / "ringfit.sh"),
    default_input_column="rfdiffusion3_path",
    build_manifest_fn=build_ringfit_manifest,
    add_run_args_fn=add_run_ringfit_args,
    collect_fn=collect_ringfit,
)
