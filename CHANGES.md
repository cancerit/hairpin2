# Changelog

## 2.0.0
  - introduced DVF flag, allowing for detection of variants arising from PCR stutter
  - fixed bug in ALF flag where variants with only 1 supporting read were not being tested, raising sensitivity
  - modified interface to streamline and promote use of config
  - changed internal architecture to a more modular approach in anticipation of further extension

## 1.2.0
  - modified calculation of statistical dispersion to match original paper, i.e. use true MAD. This gives a small increase in sensitivity

## 1.1.0

  - extended name mapping flag capabilites/syntax (see helptext and README)
