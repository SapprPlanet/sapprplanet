# sapprplanet
S-approximation method for physical fields of planets

## Description

**sapprplanet** provides several functions and example scripts for physical fields of planets.
The library allows you to make samples of satellite measurements from large data sets (currently, we use magnetic field measurements by MGS and MAVEN satellite missions, written as "*.sts" files) on arbitrary grids according to various criteria, as well as provide approximations and analytical extensions of the magnetic field (including bringing it to the same height and continuing towards sources, to the surface of the planet) with a specially developed version of the combined method for this project local and regional S-approximations, which is one of the options for implementing the method of linear integral representations.

#### Contribute
This code is still under development and benchmarking. If you find any bugs or errors in the code, please report them in GitHub or email us.

## Methods
`malkinsphere` implements method used to divide the spherical surface into equal-area cells was the spherical, rectangular equal-area grid (SREAG).
The description of the method is in the following article: 
[Malkin Z. A New Equal-area Isolatitudinal Grid on a Spherical Surface // The Astronomical J. 2019. V. 158. Iss. 4. Art. No. 158. DOI: 10.3847/1538-3881/ab3a44.](https://doi.org/10.3847/1538-3881/ab3a44)

## Model
`model` folder contains global models in spherical harmonics expansions.

`clm_Br_Maven_Plus.txt` - a model of the magnetic field of Mars, based on MAVEN data and supplemented by MGS data in regions above 75° N and below 75° S, presented as a spherical harmonic expansion.

`spectr.ipynb` - spectr of the model "Br_Maven_Plus"

## Example scripts
`malkinsphere_example` example of the malkinsphere method.

## How to install and run
If you would like to modify the source code, download the sapprplanet repository and install using pip (or pip3 depending on your installation).
```bash
    git clone https://github.com/SapprPlanet/sapprplanet.git
    cd SapprPlanet/
    pip install .
```
Alternatively, you can install sapprplanet via pip
```bash
   pip install sapprplanet
```

## To run the example scripts
```bash
    cd examples
    python malkinsphere_example.py
```

## Authors
[Alexey Batov](https://www.ipu.ru/node/82) (batov@ipu.ru),
[Anton Salnikov](https://www.ipu.ru/staff/salnikov) (salnikov@ipu.ru)

## Acknowledgments
The development of this library was supported by Russian Science Foundation grant number 23-27-00392 (https://rscf.ru/project/23-27-00392/).
