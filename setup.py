from setuptools import setup, find_packages

setup(name = 'coevo',
        version = '0.2.0',
        description = 'library for comparing coevolution identifying methods',
        url = 'in progress',
        author = 'Aram Avila-Herrera',
        author_email = 'Aram.Avila-Herrera@ucsf.edu',
        license = 'GPL',
        packages = find_packages(),
        install_requires = [ 'pandas', 'biopython' ],
        zip_safe = False
        )
