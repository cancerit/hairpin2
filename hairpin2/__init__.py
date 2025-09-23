import importlib.metadata as importlib_metadata

def _set_version() -> str:
    """Set the package version from the project metadata in pyproject.toml."""
    from warnings import warn

    fallback_version = "0.0.0"

    try:
        # __package__ allows for the case where __name__ is "__main__"
        version = importlib_metadata.version(__package__ or __name__)
    except importlib_metadata.PackageNotFoundError:
        version = fallback_version

    if version == fallback_version:
        msg = (
            f"Package version will be {fallback_version} because Python could not find "
            f"package {__package__ or __name__} in project metadata. Either the "
            "version was not set in pyproject.toml or the package was not installed. "
            "If developing code, please install the package in editable "
            "mode with `poetry install` or `pip install -e .`"
        )
        warn(msg)
    return version


__version__ = _set_version()
