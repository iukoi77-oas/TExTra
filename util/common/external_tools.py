"""Resolve vendored external tool directories."""

import glob
import os


def _parent_dirs(path):
    """Yield a path and each parent directory up to the filesystem root."""
    current = os.path.abspath(path)
    if os.path.isfile(current):
        current = os.path.dirname(current)
    while True:
        yield current
        parent = os.path.dirname(current)
        if parent == current:
            break
        current = parent


def _append_unique(items, value):
    value = os.path.abspath(value)
    if value not in items:
        items.append(value)


def external_dir_candidates(anchor_file=None):
    """Return candidate ``util/external`` directories in preferred order."""
    candidates = []

    explicit_external = os.environ.get("TEXTRA_EXTERNAL_DIR")
    if explicit_external:
        _append_unique(candidates, explicit_external)

    textra_home = os.environ.get("TEXTRA_HOME")
    if textra_home:
        _append_unique(candidates, os.path.join(textra_home, "util", "external"))

    for root in _parent_dirs(os.getcwd()):
        _append_unique(candidates, os.path.join(root, "util", "external"))

    if anchor_file:
        for root in _parent_dirs(anchor_file):
            _append_unique(candidates, os.path.join(root, "util", "external"))

    return candidates


def resolve_external_dir(anchor_file=None):
    """Return the first existing external tool directory, or the first candidate."""
    candidates = external_dir_candidates(anchor_file)
    for path in candidates:
        if os.path.isdir(path):
            return path
    return candidates[0] if candidates else os.path.abspath(os.path.join("util", "external"))


def resolve_external_file(subdir_glob, filename, anchor_file=None):
    """Return the first matching external file and all searched external roots."""
    candidates = external_dir_candidates(anchor_file)
    for external_dir in candidates:
        matches = sorted(glob.glob(os.path.join(external_dir, subdir_glob, filename)))
        for path in matches:
            if os.path.isfile(path):
                return os.path.abspath(path), candidates
    return None, candidates
