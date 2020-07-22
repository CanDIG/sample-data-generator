import setuptools

with open('README.md', 'r') as file:
    long_description = file.read()

with open('requirements.txt', 'r') as file:
    install_requires = file.readlines()

    for i in range(len(install_requires)):
        install_requires[i] = install_requires[i].strip()

setuptools.setup(
    name='sample_data_generator',
    author='Zoltan Bozoky',
    description='A script that generates 3 x 300 clinical metadata, pipeline metadata and corresponding variants for '
                'testing purposes',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/CanDIG/sample-data-generator',
    package_dir={'': 'src'},
    packages=setuptools.find_packages(where='src'),
    python_requires='>=3.7',
    install_requires=install_requires,
    entry_points={
        'console_scripts': [
            'generate_sample_data=sample_data_generator.generate_vcf:main',
        ],
    },
)
