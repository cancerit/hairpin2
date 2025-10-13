class ConstraintError(ValueError):
    pass


def bound[T: float | int](val: T, min: T | None = None, max: T | None = None) -> T:
    min_check = min <= val if min is not None else True
    max_check = val <= max if max is not None else True
    if not (min_check and max_check):
        raise ConstraintError(
            f"input value {val} is outside of bound range min: {min} - max: {max}"
        )
    return val
