from setuptools import setup

exec(compile(open('SemiBin/semibin_version.py').read(),
             'SemiBin/semibin_version.py', 'exec'))

long_description = open('README.md', encoding='utf-8').read()

setup(name='SemiBin',
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
        'Programming Language :: Python :: 3.9',
      ],
      url='https://github.com/BigDataBiology/SemiBin',
      author='Shaojun Pan',
      author_email='shaojun1997777@gmail.com',
      license='MIT',
      packages = ['SemiBin'],
      setup_requires=['numpy'],
      install_requires=[
          'Biopython',
          'tqdm',
          'pyyaml',
          'atomicwrites',
          'torch==1.8',
          'python-igraph',
          'pandas==1.2.4',
          'scikit-learn',
          'requests',
      ],
      package_data={
          'SemiBin': ['*.hmm','*.pl','*.h5']},
      zip_safe=False,
      entry_points={
            'console_scripts': ['SemiBin=SemiBin.main:main'],
      }
)

