from hairpin2 import abstractfilters as haf
from collections.abc import Collection
from typing import Any
# pyright: reportExplicitAny=false
# pyright: reportAny=false

type AnyResultCollection = Collection[haf.FilterResult[Any]]
