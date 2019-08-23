import logging

logging.basicConfig(format = u'%(levelname)-8s [%(asctime)s] %(message)s', level = logging.DEBUG)
try:
    from setuptools import setup, find_packages
except:
    logging.error('setuptools is not installed or outdated!\n\n'
                 'You can install or update setuptools using\n'
                 'pip install --upgrade setuptools (if you have pip)\n'
                 'or\n'
                 'sudo apt-get install python-setuptools (on Ubuntu)\n')
    exit(-1)

setup(name='orf-search',
      version='0.4',
      author='Tatiana Dvorkina',
      author_email='tanunia@gmail.com',
      url='https://github.com/ablab/orf-search',
      license='GPLv2',
      packages=find_packages(),
      install_requires=[
        'joblib',
        'edlib',
        'argparse',
        'pyyaml',
        'biopython'
      ],
      include_package_data = True
      )
