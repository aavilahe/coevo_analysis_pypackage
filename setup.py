from setuptools import setup

setup(name = 'coevo',
        version = '0.1',
        description = 'library for comparing coevolution identifying methods',
        url = 'in progress',
        author = 'Aram Avila-Herrera',
        author_email = 'Aram.Avila-Herrera@ucsf.edu',
        license = 'GPL',
        packages = [ 'coevo' ],
        install_requires = [ 'pandas', 'biopython' ],
        zip_safe = False
        )
