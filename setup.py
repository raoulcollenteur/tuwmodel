from setuptools import setup, find_packages

# Get the version.
version = {}
with open("tuwmodel/version.py") as fp:
    exec(fp.read(), version)

setup(
    name='tuwmodel',
    version=version['__version__'],
    platforms='Windows, Mac OS-X',
    url='http://github.com/raoulcollenteur/tuwmodel',
    license='MIT License',
    author='Raoul Collenteur',
    author_email='raoulcollenteur@gmail.com',
    description='python implementation of the rainfall runoff model TUWmodel '
                'developed by TU Wien (Parajka et al, 2007)',
    packages=find_packages(exclude=[]),
    package_data={"tuwmodel": ["hbvmodel.so"], },
)
