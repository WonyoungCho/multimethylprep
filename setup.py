# Lib
from setuptools import setup, find_packages
exec(open('multimethylprep/version.py').read())

test_requirements = [
    'methylcheck', # 'git+https://github.com/FoxoTech/methylcheck.git@feature/v0.7.7#egg=methylcheck',
    'pytest',
    'pytest_mock',
    'matplotlib',
    'scikit-learn', # openpyxl uses this, and forcing it to install the best version, not sklearn 0.0
    'openpyxl',
    'coverage'
]

setup(
    name='multimethylprep',
    version=__version__,
    description='Python-based Illumina methylation array preprocessing software',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    project_urls = {
        "Documentation": "https://life-epigenetics-methylprep.readthedocs-hosted.com/en/latest/",
        "Source": "https://github.com/WonyoungCho/methylprep/",
    },
    python_requires='>=3.6',
    keywords='methylation dna data processing epigenetics illumina',
    url='https://github.com/WonyoungCho/multimethylprep',
    license='MIT',
    author='Wonyoung Cho',
    author_email='bourbaki10@gmail.com',
    packages=find_packages(),
    install_requires=[
        'pyparsing > 3.0',
        'numpy',
        'pandas >=1.3.0',
        'scipy',
        'statsmodels',
        'parmap',
        'tqdm',
        'bs4',
        'lxml',
        'requests',
    ],
    extras_require={
        'dev': test_requirements
    },
    setup_requires=['pytest-runner'],
    tests_require= test_requirements,
    entry_points='''
        [console_scripts]
        methylprep-cli=methylprep.cli:cli_app
    ''',
)
