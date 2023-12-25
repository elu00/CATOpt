To install dependencies:
```
sudo apt install libsuitesparse-dev libeigen3-dev
sudo apt install gfortran libopenblas-dev liblapack-dev liblapacke-dev libmetis-dev 
sudo apt install intel-mkl
```
Cloning and build:
```
git clone --recursive git@github.com:elu00/CATOpt.git
cd CATOpt
mkdir build
cd build
cmake ..
make -j12
```
