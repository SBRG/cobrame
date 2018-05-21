from setuptools import setup, find_packages


setup(name="cobrame",
      version="0.0.9",
      author="Colton Lloyd and Ali Ebrahim",
      url="https://github.com/SBRG/cobrame",
      install_requires=["sympy", "six", "Biopython", "cobra<=0.5.11", "pandas",
                        "scipy", "numpy", "setuptools", "jsonschema"],
      package_data={'': ['io/JSONSCHEMA']},
      packages=find_packages()
      )
