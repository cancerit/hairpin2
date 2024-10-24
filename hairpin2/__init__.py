# hairpin2
#
# Copyright (C) 2024 Genome Research Ltd.
#
# Author: Alex Byrne <ab63@sanger.ac.uk>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.

# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


def _set_version() -> str:  # noqa: C901
    """Set the package version from the project metadata in pyproject.toml."""
    from warnings import warn

    fallback_version = "0.0.0"
    try:
        # importlib.metadata is present in Python 3.8 and later
        import importlib.metadata as importlib_metadata
    except ImportError:
        # use the shim package importlib-metadata pre-3.8
        import importlib_metadata as importlib_metadata

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
