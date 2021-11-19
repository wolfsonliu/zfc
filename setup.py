#! /usr/bin/env python3

import os
try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

package = 'zfc'
version = '0.1.7'


def readme():
    with open(
            os.path.join(
                os.path.dirname(os.path.abspath(__file__)),
                'README.md'
            )
    ) as f:
        return f.read()

long_description=readme()

setup(
    name=package,
    version=version,
    description="The zfc software is used for analysis of counts data.",
    long_description=long_description,
    long_description_content_type='text/markdown',
    classifiers=[
        'Development Status :: 1 - Planning',
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)'
    ],
    url='https://github.com/wolfsonliu/zfc',
    author='Zhiheng Liu',
    author_email='zhiheng.liu@pku.edu.cn',
    license='GPL',
    packages=['zfc'],
    install_requires=[
        'numpy', 'scipy', 'pandas',
        'matplotlib', 'sklearn'
    ],
    scripts=[
        'bin/zfc',
        'bin/library_count_sgrna',
        'bin/library_count_sgrna_with_barcode'
    ],
    package_dir={'zfc': 'zfc'},
    include_package_data=True,
    zip_safe=False
)
