coevo
=====

A library of tools and scripts for preprocessing
sequence and structure files for coevolutionary
analyses.

----

See `SEQUENCE_FORMATTING`_ and `MAPPING_AND_DISTANCES`_ for example usage.

Install
-------

.. code-block:: bash

    pip install git+https://github.com/aavilahe/coevo_analysis_pypackage.git@dev

Scripts in ``bin/``
-------------------

.. code-block:: bash

    convert_resnums_to_columns.py
    fasta_to_phy.py
    fasta_to_psicov.py
    get_dists.py
    join_fastas.py
    make_attributes.py
    map_column_to_resnum.py
    min_dists.py
    split_faa_on_col.py

Previously, this package provided support for postprocessing
output from various coevolution programs. This functionallity
is deprecated and will be moved to a separate package in the future.

.. _SEQUENCE_FORMATTING: doc/SEQUENCE_FORMATTING.md
.. _MAPPING_AND_DISTANCES: doc/MAPPING_AND_DISTANCES.md
