#! /usr/bin/env python3

import os
try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

package = 'zfc'
version = '0.1.0'


def readme():
    with open(
            os.path.join(
                os.path.dirname(os.path.abspath(__file__)),
                'README.md'
            )
    ) as f:
        return f.read()


setup(
    name=package,
    version=version,
    description="The zfc software is used for analysis of counts data.",
    long_description=readme(),
    classifiers=[
        'Development Status :: 1 - Planning',
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)'
    ],
    url='',
    author='Zhiheng Liu',
    author_email='zhiheng.liu@pku.edu.cn',
    license='GPL',
    packages=['zfc'],
    install_requires=[
        'numpy>=1.10', 'scipy>=1.0', 'pandas>=0.16'
    ],
    scripts=['bin/zfc'],
    package_dir={'zfc': 'zfc'},
    include_package_data=True,
    zip_safe=False
)
