"""
Earth Engine deployment settings of the EII package.

In a source checkout, the ``gee`` block of the project-level ``config.yaml`` sets the Earth
Engine project, the root of the v1 asset paths and the default dataset version. Installed
wheels do not ship ``config.yaml`` and use the built-in defaults below.
"""

from pathlib import Path

import yaml  # type: ignore[import-untyped]

# Determine project root
# Assuming this file is in src/eii/, the root is 2 levels up
PROJECT_ROOT = Path(__file__).parent.parent.parent
CONFIG_PATH = PROJECT_ROOT / "config.yaml"


def load_config():
    """Load configuration from config.yaml or return defaults."""
    if CONFIG_PATH.exists():
        with open(CONFIG_PATH) as f:
            return yaml.safe_load(f)
    return {}


_config = load_config()
_gee_config = _config.get("gee", {})

# Public Constants
GEE_PROJECT = _gee_config.get("project", "landler-open-data")
EII_ASSET_ROOT = _gee_config.get("asset_root", "projects/landler-open-data/assets/eii")
VERSION = _gee_config.get("version", "v2")
