try:
    from setuptools import setup, find_packages
except:
    logger.error('setuptools is not installed or outdated!\n\n'
                 'You can install or update setuptools using\n'
                 'pip install --upgrade setuptools (if you have pip)\n'
                 'or\n'
                 'sudo apt-get install python-setuptools (on Ubuntu)\n',
                 exit_with_code=1)

setup(name='orf-search',
      version='0.4',
      author='Tatiana Dvorkina',
      author_email='tanunia@gmail.com',
      url='https://github.com/ablab/orf-search',
      license='GPLv2',
      packages=['','scripts'],
      install_requires=[
        'joblib',
        'edlib',
        'argparse',
        'pyyaml',
        'biopython'
      ],
      include_package_data = True
      )
