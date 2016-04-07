from distutils.core import setup

setup(
    name='starmod',
    version='0.0.1',
    packages=['starmod', 'starmod.src', 'starmod.src.lib'],
    scripts=['scripts/lineProf.py'],
    url='http://github.com/tribeiro/starmod',
    license='GPL v2',
    author='Tiago Ribeiro',
    author_email='tribeiro@ufs.br',
    description='Star modeling tools.'
)
