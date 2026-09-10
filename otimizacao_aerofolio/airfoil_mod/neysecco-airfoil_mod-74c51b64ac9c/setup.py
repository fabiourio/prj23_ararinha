from setuptools import setup, find_packages

setup(
    name='airfoil_mod',
    version='0.1',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'scipy',
        'matplotlib'
    ],
    author='Ney Rafael Sêcco',
    author_email='ney.secco@gp.ita.br',
    description='Module that generates and manipulates airfoil sections',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    #url='https://github.com/yourusername/airfoil_mod',  # Project URL
    #classifiers=[
    #    'Programming Language :: Python :: 3',
    #    'License :: OSI Approved :: MIT License',
    #    'Operating System :: OS Independent',
    #],
    python_requires='>=3.6',
)
