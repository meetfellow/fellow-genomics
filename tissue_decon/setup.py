from setuptools import setup, find_packages

def read_requirements():
    with open('requirements.txt') as req:
        content = req.read()
        requirements = content.split('\n')
    requirements = [r for r in requirements if r != '']
    return requirements


setup(
    name='tissue-decon',
    version='0.1.0',
    description='Tools for deconvoluting tissue methylation data',
    author='James Hart',
    author_email='jhart@meetfellow.com',
    packages=find_packages(exclude=('tests', 'docs')),
    install_requires=read_requirements(),
    entry_points={
    'console_scripts': [
        'tissue-methyl-decon=tissue_decon.tissue_methyl_decon:main'
    ]
}
)