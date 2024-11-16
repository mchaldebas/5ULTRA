from setuptools import setup, find_packages
import io
import os

# Read the contents of your README file
current_directory = os.path.abspath(os.path.dirname(__file__))
with io.open(os.path.join(current_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(name='5ULTRA',
        description='5ULTRA: A Pipeline for 5UTR Variant Annotation and Scoring',
        long_description=long_description,
        long_description_content_type='text/markdown',
        author='Matthieu Chaldebas',
        author_email='mchaldebas@rockefeller.edu',
        license='GPLv3',
        url='https://github.com/mchaldebas/5ULTRA',
        packages=['5ULTRA'],
        install_requires=['pysam>=0.22.1',
                            'pandas>=1.4.4',
                            'numpy>=1.21.5',
                            'joblib>=1.1.0',
                            'scikit-learn>=1.0.2',
        ],
        package_data={'5ULTRA': ['data/*',
                                'data/5UTR.hg38.phastCons100way/*',
                                'data/5UTR.hg38.phyloP100way/*',
                                ],
        },
        entry_points={'console_scripts': [
                        '5ULTRA=5ULTRA.__main__:main',
                        ],
        },
        classifiers=['Programming Language :: Python :: 3',
                    'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
                    'Operating System :: OS Independent',
                    ],
        python_requires='>=3.6',
        test_suite='tests',
        )
