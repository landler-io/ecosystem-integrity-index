"""
Configuration utility to load key=value config files.
"""

import contextlib
from pathlib import Path
from typing import Any


def load_config(path: Path) -> dict[str, Any]:
    """
    Parses a simple configuration file with key=value syntax.
    Supports comments starting with #.
    Infers types: int, float, boolean, None.
    Strings can be quoted or unquoted.
    """
    config = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            if "=" in line:
                key, val = line.split("=", 1)
                key = key.strip()
                val = val.strip()

                # Handle inline comments if not quoted
                # Simple check: if quote is at start, look for quote at end.
                # If content is inside quotes, don't strip #.
                is_quoted = (val.startswith('"') and val.endswith('"')) or (
                    val.startswith("'") and val.endswith("'")
                )

                if not is_quoted and "#" in val:
                    val = val.split("#", 1)[0].strip()

                # Strip quotes
                if is_quoted:
                    val = val[1:-1]
                else:
                    # Type inference
                    val_lower = val.lower()
                    if val_lower == "true":
                        val = True
                    elif val_lower == "false":
                        val = False
                    elif val_lower in ("none", "null"):
                        val = None
                    else:
                        with contextlib.suppress(ValueError):
                            val = float(val) if "." in val else int(val)

                config[key] = val
    return config
