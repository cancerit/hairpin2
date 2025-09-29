# Changes

2.0.0 -> 3.0.0

- Changes to {py:meth}`~hairpin2.sci_funcs.tag_stutter_duplicates()`:
    - Now selects highest quality read as representative of a duplicate group. This is important and influential on flag calls since it changes the overlap between low quality tags and duplicate tags on reads.
