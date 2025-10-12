.. hairpin2 documentation master file, created by
   sphinx-quickstart on Sat Sep 27 11:10:01 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

hairpin2 documentation
======================

``hairpin2`` – **read-aware artefactual variant flagging** - GitHub_

hairpin2 is designed to flag variants that are likely artefactual via a series of tests performed upon the read data associated with each variant. Initially, it was concieved to flag possible cruciform artefacts for LCM sequence data, but the concept has been extended and can detect a variety of potentially spurious variants (including indels). The tool operates on a VCF file containing one or more samples, and alignment files for all samples to be tested.

Read :ref:`quickstart-target` to get started!

.. _GitHub: https://github.com/cancerit/hairpin2

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   quickstart
   guide
   tests
   hairpin2
   dev
   changes
