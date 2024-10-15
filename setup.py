from setuptools import setup

setup(
    name='aerospike',
    version='0.1',
    packages=['src'],
    install_requires=[
            "numpy",
            "matplotlib",
            "scipy",
            "CoolProp",
            "pint",
            "icecream",
            "logging",
      ],
)