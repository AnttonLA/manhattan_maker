#!/usr/bin/env python

from setuptools import setup


VERSION = "1.7.4"


def setup_package():
    setup(
        name="manhattan_maker",
        version=VERSION,
        description="Creation of beautiful Manhattan plots",
        author="Louis-Philippe Lemieux Perreault",
        author_email="louis-philippe.lemieux.perreault@statgen.org",
        maintainer="Antton Lamarca",
        maintainer_email="antton.lamarca@med.lu.se",
        url="https://github.com/AnttonLA/manhattan_maker",
        license="CC BY-NC 4.0",
        entry_points={
            "console_scripts": [
                "manhattan_maker=manhattan_maker:safe_main",
            ],
        },
        py_modules=["manhattan_maker"],
        install_requires=["matplotlib >=1.3.1", "numpy >= 1.8.0",
                          "pandas >= 0.17.0"],
        classifiers=[
            "Operating System :: Linux",
            "Programming Language :: Python",
            "Programming Language :: Python :: 2.7",
        ],
    )

    return


if __name__ == "__main__":
    setup_package()
