from setuptools import setup, find_packages

setup(
    name='dimer_calculations',
    version='0.1.0',
    packages=find_packages(),
    install_requires=[
        'numpy>=1.18.0',
        'scipy>=1.4.0',
        # Add other dependencies here
    ],
    entry_points={
        'console_scripts': [
            # Define any command-line scripts here, for example:
            # 'dimer_calculations=dimer_calculations:main',
        ],
    },
    author='Your Name',
    author_email='your.email@example.com',
    description='A short description of your project',
    url='https://github.com/ewolpert1/dimer_calculations',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.11',
)
