import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="ga_variable_length-IgnacioGarridoBotella", # Replace with your own username
    version="0.0.2",
    author="Ignacio Garrido Botella",
    author_email="ignaciogabo95@gmail.com",
    description="Proyect to launch a genetic algorithm.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/IgnacioGarrido/ga_variable_length.git",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ]
    )