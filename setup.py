try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

setup(name="cobrame",
      version="0.0.5",
      author="Ali Ebrahim and Colton Lloyd",
      author_email="minime_dev@googlegroups.com",
      url="https://github.com/SBRG/cobrame",
      install_requires=["sympy", "six", "Biopython", "cobra", "pandas",
                        "scipy", "numpy"],
      packages=["cobrame"]
      )
