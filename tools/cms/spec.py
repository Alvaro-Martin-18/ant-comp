from pathlib import Path

from prosapia.core import Tool

from .collect_cms import collect_cms
from .run_cms import NO_DEFAULT_COLUMN, add_run_cms_args, build_cms_manifest

TOOL = Tool(
    name="cms",
    action="update",
    description="Contact molecular surface + shape complementarity (GPU, cms-cuda).",
    default_script=str(Path(__file__).parent / "cms.sh"),
    # Sentinel, as in usalign/chainsel: no structure column is a defensible default,
    # so the manifest builder raises unless -i/--input-column is given.
    default_input_column=NO_DEFAULT_COLUMN,
    build_manifest_fn=build_cms_manifest,
    add_run_args_fn=add_run_cms_args,
    collect_fn=collect_cms,
)
