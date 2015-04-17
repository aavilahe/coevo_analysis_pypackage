from setuptools import setup, find_packages

setup(name = 'coevo',
        version = '0.2.0',
        description = 'library for comparing coevolution identifying methods',
        url = 'in progress',
        author = 'Aram Avila-Herrera',
        author_email = 'Aram.Avila-Herrera@ucsf.edu',
        license = 'GPL',
        packages = find_packages(),
        install_requires = ['pandas', 'biopython'],
        scripts = [
                   'bin/fasta_to_phy.py',
                   'bin/fasta_to_psicov.py',
                   'bin/get_dists.py',
                   'bin/join_fastas.py',
                   'bin/make_attributes.py',
                   'bin/map_column_to_resnum.py',
                   'bin/min_dists.py',
                   'bin/split_faa_on_col.py'
                   ],
        zip_safe = False
        )
