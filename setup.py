import setuptools
import glob
import os
from sccore.__init__ import VERSION

with open("README.md", "r") as fh:
    long_description = fh.read()


def parse_requirements(filename):
    with open(filename, "r") as file:
        lines = file.readlines()
    requirements = [line.strip() for line in lines if line.strip() and not line.startswith("#")]
    return requirements


entrys = []
for file in glob.glob("sccore/cli/*.py"):
    name = os.path.basename(file)[:-3]
    cli = name.replace("_", "-")
    entrys.append(f"{cli}=sccore.cli.{name}:main")
entry_dict = {
    "console_scripts": entrys,
}

setuptools.setup(
    name="sccore",
    version=VERSION,
    author="zhouyiqi",
    author_email="zhouyiqi@singleronbio.com",
    description="Single Cell core functions",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/singleron-RD/sccore",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.9",
    include_package_data=True,
    entry_points=entry_dict,
    install_requires=parse_requirements("requirements.txt"),
)
