"""Read and update the global TExTra run configuration."""

import json
import os


GLOBAL_CONFIG_NAME = "TExTra.config.json"


def global_config_path(out_dir):
    """Return the global run config path under an output directory."""
    return os.path.join(out_dir, GLOBAL_CONFIG_NAME)


def read_global_config(out_dir):
    """Read ``TExTra.config.json`` from an output directory, or return an empty config."""
    path = global_config_path(out_dir)
    if not os.path.isfile(path):
        return {}, path
    with open(path, "r", encoding="utf-8") as handle:
        payload = json.load(handle)
    if not isinstance(payload, dict):
        raise RuntimeError(f"Invalid TExTra config: {path} must contain a JSON object.")
    return payload, path


def update_global_config(out_dir, module_name, module_payload):
    """Merge one module payload into the global run config."""
    config, path = read_global_config(out_dir)
    modules = config.get("modules", {})
    if not isinstance(modules, dict):
        modules = {}
    abs_out_dir = os.path.abspath(out_dir)
    config.update(
        {
            "schema_version": 1,
            "output_root": abs_out_dir,
            "modules": modules,
        }
    )
    modules[module_name] = dict(module_payload)
    os.makedirs(out_dir, exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(config, handle, indent=2, sort_keys=True)
    return path


def extract_module_config(config, module_name):
    """Return one module config from a loaded global config."""
    modules = config.get("modules", {})
    if not isinstance(modules, dict):
        return {}
    module_config = modules.get(module_name, {})
    return module_config if isinstance(module_config, dict) else {}
