from setuptools import setup, find_packages

setup(
  name='pyexocross',
  version='0.0.1',
  description='A Python wrapper for the Exomol database and the ExoCross code.',
  long_description='...',  # Optional: Longer description
  long_description_content_type="text/markdown",  # Optional: If using markdown
  url='https://github.com/DGonzalezPicos/pyexocross',
  author='Dario Gonzalez Picos',
  author_email='picos@strw.leidenuniv.nl',
  license='MIT',
  install_requires=[],  # Optional: List of external dependencies
  packages=find_packages(),
  include_package_data=True  # Optional: If you have data files in your package
)
