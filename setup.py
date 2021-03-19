from setuptools import setup

exec(compile(open('s3n2bin/s3n2bin_version.py').read(),
             's3n2bin/s3n2bin_version.py', 'exec'))

long_description = open('README.md', encoding='utf-8').read()

setup(name='S3N2Bin',
      version=__version__,
      description='Metagenomic binning with semi-supervised siamese neural network',
      long_description = long_description,
      long_description_content_type = 'text/markdown',
      classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Artificial Intelligence',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
      ],
      url='https://github.com/BigDataBiology/S3N2Bin',
      author='Shaojun Pan',
      author_email='shaojun1997777@gmail.com',
      license='MIT',
      packages = ['s3n2bin'],
      install_requires=[
          'numpy',
          'Biopython',
          'tqdm',
          'pyyaml',
          'atomicwrites',
          'torch',
          'python-igraph'
      ],
      package_data={
          's3n2bin': ['*.hmm','*.pl']},
      zip_safe=False,
      entry_points={
            'console_scripts': ['S3N2Bin=s3n2bin.main:main'],
      }
)

