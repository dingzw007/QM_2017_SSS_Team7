[![codecov](https://codecov.io/gh/MolSSI-SSS/QM_2017_SSS_Team7/branch/master/graph/badge.svg)](https://codecov.io/gh/MolSSI-SSS/QM_2017_SSS_Team7)
[![Build Status](https://travis-ci.org/MolSSI-SSS/QM_2017_SSS_Team7.svg?branch=master)](https://travis-ci.org/MolSSI-SSS/QM_2017_SSS_Team7)

# QM_2017_SSS_Team7

## Installation
Download the repository and run `setup.py` as follows:
```
git clone https://github.com/MolSSI-SSS/QM_2017_SSS_Team7
cd QM_2017_SSS_Team7
python setup.py install
```
Alternatively you can set up the environment by following the instructions [here](https://molssi-sss.github.io/Logistics_SSS_2017/Setup.html). After the environment is configures make sure to activate your python environment:

```python
source activate sss
```

## Usage
Run following to see options:
```
python scf.py --help
```

### Compiling the JK algorithm
Go to `JKbinding` folder and run:
```
mkdir build
cmake ../ -DCMAKE_CXX_FLAGS="-Ofast"
cd build
make
```
Then copy the Python code to `build` directory and run it:
```
cp ../JKbinding.py .
python JKbinding.py
```

## Team Members
- Srinivas Mushnoori
- Sangeeta Sur
- Zhiwei Ding
- Kutay B. Sezginel
