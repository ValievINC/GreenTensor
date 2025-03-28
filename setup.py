from setuptools import setup, find_packages

setup(
    name="GreenTensor_Extended",            
    version="0.1",                
    author="Denisov Dmitry, Valiev Bulat",           
    description="library for the analysis of electromagnetic field scattering on heterogeneous spherical structures.", 
    long_description=open("README.rst").read(),
    long_description_content_type="text/x-rst",
    packages=find_packages(),
    install_requires=[
        "matplotlib",
        "scipy",
        "numpy"
    ],
    python_requires=">=3.12",
    data_files=[("", ["example.py"])],
)