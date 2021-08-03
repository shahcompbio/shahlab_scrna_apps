from setuptools import setup, find_packages

setup(
    name='gex_analysis',
    version='0.0.1',
    packages=['gex_analysis'],
    py_modules=['gex_analysis.gex_analysis'],
    long_description=open('README.md').read(),
    entry_points={
        'console_scripts': [
            'gex_analysis = gex_analysis.gex_analysis:cli',
        ],
    },
)