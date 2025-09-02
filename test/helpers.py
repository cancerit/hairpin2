from dataclasses import fields, is_dataclass
from typing import Any
from hairpin2.structures import ReadView

def comp_ReadView(s1: ReadView, s2: ReadView):
    return all(a == b for a, b in zip(s1.all, s2.all)) and len(s1.all) == len(s2.all) and len(s1) == len(s2)


def unsafe_construct(cls: type, **data: Any):
    """
    Build a dataclass instance without running Pydantic/Dataclass validation,
    bypassing `frozen=True` by using object.__setattr__.
    """
    if not is_dataclass(cls):
        raise TypeError("unsafe_construct expects a dataclass type")

    inst = object.__new__(cls)  # no __init__, no validation
    valid_names = {f.name for f in fields(cls)}

    # set declared dataclass fields, bypass frozen
    for name, value in data.items():
        if name in valid_names:
            object.__setattr__(inst, name, value)

    return inst


