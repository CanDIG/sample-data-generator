# Randomly generated 3 x 300 patient data

For testing purposes 3 x 300 clinical metadata, pipeline metadata and corresponding variants have been generated.

Check releases for data files!

## Building
To build, you must have `wheel` installed. If you do not have it installed, run the following:

`pip install wheel`

Once you have installed `wheel`, run the following at the root level of this repository to build:

`python setup.py sdist bdist_wheel`

After this command has successfully run, a directory called `dist` should exist. Within this directory 
there should be a `.whl` file and a `.tar.gz` file.

## Installing
To install, you have two options available:
1. Install from a built distribution
2. Install from source code

### Installing From a Built Distribution
To install a built distribution, you must have `wheel` installed. If you do not have it installed, run the following:

`pip install wheel`

Once you have installed `wheel`, go into the `dist` directory and run the following to install:

`pip install sample_data_generator-0.0.0-py3-none-any.whl`

After this command has successfully run, it should be installed and ready to run.

### Installing From Source Code
To install, run the following at the root level of this repository:

`pip install .`

After this command has successfully run, it should be installed and ready to run.

## Running
```
Generate random vcf files based on resampled dbSNP data

Usage:                                                                            
  generate_sample_data <number_of_files> <output_folder> <filename_start>

Options:
  -h --help          Show this screen                                                                                     
  <number_of_files>  Number of samples to generate                                                                        
  <output_folder>    Path information where the generated files will be put                                               
  <filename_start>   All filename will start with this phrase and continued by                                                               
                     a number.                                                                                                                                                                                                                  

Examples:                                                                                                                 
  generate_sample_data 5000 tmp/variants/ patient_
```

## Development
Create a virtual environment:

`python3 -m venv venv`

Activate the virtual environment:

`source venv/bin/activate`

Go to the root level of this repository and install it in editable mode:

`pip install -e .`