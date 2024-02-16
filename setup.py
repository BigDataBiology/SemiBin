from setuptools import setup

exec(compile(open('SemiBin/semibin_version.py').read(),
             'SemiBin/semibin_version.py', 'exec'))

long_description = open('README.md', encoding='utf-8').read()

setup(name='SemiBin',
      version=__version__,
      python_requires='>=3.7',
      description='Metagenomic binning with siamese neural networks',
      long_description = long_description,
      long_description_content_type = 'text/markdown',
      classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Artificial Intelligence',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
      ],
      url='https://github.com/BigDataBiology/SemiBin',
      author='Shaojun Pan',
      author_email='shaojun1997777@gmail.com',
      license='MIT',
      packages = ['SemiBin'],
      install_requires=[
          'numpy',
          'tqdm',
          'pyyaml',
          'torch',
          'python-igraph',
          'pandas',
          'scikit-learn',
          'requests',
          'numexpr',
      ],
      package_data={
          'SemiBin': ['*.hmm', 'models/*.pt']},
      zip_safe=False,
      entry_points={
            'console_scripts': ['SemiBin=SemiBin.main:main_no_version',
                                'SemiBin2=SemiBin.main:main2',
                                'SemiBin1=SemiBin.main:main1',
                                ],
      }
)

