try:
    from dune.packagemetadata import metaData, Description
except ImportError:
    from packagemetadata import metaData, Description
from setuptools import find_namespace_packages
from skbuild import setup

description = Description('dune.module')
version = description.versionstring.replace('-git', '')

buildVersion = '.0.dev20201212'
duneVersion = '2.8' + buildVersion

metadata = metaData(duneVersion)[1]
metadata['version'] = version + buildVersion
metadata['long_description'] = metadata['long_description'].replace('doc/logo/dumux_logo_hires_whitebg.png', 'https://dumux.org/images/logo.svg')
metadata['packages'] = find_namespace_packages(where='python', include=['dumux.*'])

setup(**metadata)
