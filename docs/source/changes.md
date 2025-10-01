# Change Log

## 2.0.0 -> 3.0.0

  - Changes to {py:meth}`~hairpin2.sci_funcs.tag_stutter_duplicates()`:
    - Now selects highest quality read as representative of a duplicate group. This is important and influential on flag calls since it changes the overlap between low quality tags and duplicate tags on reads.

## 1.2.0 -> 2.0.0
  - introduced DVF flag, allowing for detection of variants arising from PCR stutter
  - fixed bug in ALF flag where variants with only 1 supporting read were not being tested, raising sensitivity
  - modified interface to streamline and promote use of config
  - changed internal architecture to a more modular approach in anticipation of further extension

## 1.1.0 -> 1.2.0
  - modified calculation of statistical dispersion to match original paper, i.e. use true MAD. This gives a small increase in sensitivity

## 1.0.0 -> 1.1.0

  - extended name mapping flag capabilites/syntax (see helptext and README)