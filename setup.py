try:
    from dune.packagemetadata import metaData
except ImportError:
    from packagemetadata import metaData
from setuptools import find_namespace_packages
from skbuild import setup

# When building a new package, update the version numbers below and run:
# > python setup.py sdist
# > python -m twine upload dist/*

dumuxVersion = '3.4.0.dev20210114'
duneVersion  = '2.8.0.dev20201218'

metadata = metaData(duneVersion)[1]
metadata['version'] = dumuxVersion
metadata['long_description'] = metadata['long_description'].replace(
    'doc/logo/dumux_logo_hires_whitebg.png',
    'https://dumux.org/images/logo.svg'
)
metadata['packages'] = find_namespace_packages(where='python', include=['dumux.*'])

# auto-generate pyproject.toml with duneVersion when building sdist
from skbuild.command.sdist import sdist
class mysdist(sdist):
    def run(self):
        requires = ['setuptools', 'wheel', 'scikit-build', 'cmake', 'ninja', 'requests']
        requires += metadata['install_requires']
        with open('pyproject.toml', 'w') as f:
            f.write("[build-system]\n")
            f.write("requires = ['"+"', '".join(requires)+"']\n")
            f.write("build-backend = 'setuptools.build_meta'\n")
        sdist.run(self)
metadata['cmdclass'] = {'sdist': mysdist}

setup(**metadata)
