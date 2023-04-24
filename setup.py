from setuptools import setup, find_packages

setup(
    name="atomate-tsai",
    author="Jeng-Yuan Tsai",
    version="0.1.0",
    author_email="tsaie79@gmail.com",
    summary="Atomate fork for Tsai",
    packages=find_packages(),
    install_requires=[
        "smoqe@git+https://github.com/dangunter/smoqe.git@master#egg=smoqe"
        "pymatgen_db@git+https://github.com/materialsproject/pymatgen-db.git@c3271276c2ef26dc98ccc86634405a04cd677395"
        "#egg=pymatgen_db",
        "pymatgen_diffusion@git+https://github.com/materialsvirtuallab/pymatgen-analysis-diffusion.git"
        "@51d84ea1fd034941baa22d4b1b8610a4cc6fb801#egg=pymatgen_diffusion",
        "ase@git+ssh://git@github.com/rosswhitfield/ase.git@07de35654601ddbb2b23a4e7df7091696b0af108#egg=ase",
        "pytopomat@git+https://github.com/tsaie79/pytopomat.git@38e5856eec61800a52b44ea0dfd2c99de311f97f#egg=pytopomat",
        "pycdt@git+https://github.com/tsaie79/pycdt.git@master#egg=pycdt",
        "pymatgen@git+https://github.com/tsaie79/pymatgen.git@master#egg=pymatgen",
        "atomate@git+https://github.com/tsaie79/atomate.git@master#egg=atomate",
        "FireWorks==1.9.6",
        "custodian==2020.4.27",
        "phonopy==2.12.0",
        "pymongo==3.11.0",
    ],
)
