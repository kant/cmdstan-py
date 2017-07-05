from setuptools import setup, find_packages

__version__ = "0.1.1"

setup(
    name="cmdstan-py",
    version=__version__,
    packages=find_packages(),
    install_requires=[
        'six',
        'numpy',
        'pandas'
    ],
    # command
    scripts=['BayesGLM.py'],
    # metadata
    author="Linnarsson Lab",
    author_email="gioelelamanno@gmail.com",
    description="cmd-stan wrapper",
    license="MIT",
    url="https://github.com/linnarsson-lab/cmdstan-py",
)
