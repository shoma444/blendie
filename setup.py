# Created by Shoma Yamanouchi (shoma.yamanouchi@ec.gc.ca)

######## Requirements ########
# Python 3+ (3.9+ recommended) 
# numpy
# datetime
# metpy
# herbie
# grib2io
##############################

# IMPORTANT NOTICE: This package will ONLY work in North America #


import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="blendie",
    version="0.0.2",
    author="Shoma Yamanouchi",
    author_email="shoma.yamanouchi@ec.gc.ca",
    description="",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/shoma444/blendie",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)

if __name__ == "__main__":
    setuptools.setup()