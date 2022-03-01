import glob
from setuptools import setup, find_packages


with open("requirements.txt") as f:
    install_requires = [x.rstrip() for x in f]

setup(
    name="cte",
    version="0.0.1",
    description="Evaluate accuracy of Covid consensus sequence"
    packages=find_packages(exclude=["tests"]),
    package_data={'cte': ['data/*']},
    author="Martin Hunt",
    author_email="mhunt@ebi.ac.uk",
    url="https://github.com/iqbal-lab-org/covid-truth-eval",
    test_suite="nose.collector",
    tests_require=["pytest"],
    entry_points={"console_scripts": ["cte = cte.__main__:main"]},
    install_requires=install_requires,
    license="MIT",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3 :: Only",
        "License :: OSI Approved :: MIT License",
    ],
)

