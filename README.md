# How to install

Reverse-engineered from Rohan's [mpMap2](https://github.com/rohan-shah/mpMap2)
installation instructions

```
git clone https://github.com/rohan-shah/multiStateFlowReliability.git
git clone https://github.com/rohan-shah/Rcpp.git
git clone https://github.com/rohan-shah/maxFlowAlgorithms.git
cd Rcpp
mkdir release
cd release
cmake .. -DCMAKE_BUILD_TYPE=Release
make
cd ../../maxFlowAlgorithms
mkdir release
cd release
cmake .. -DCMAKE_BUILD_TYPE=Release -DRcpp_DIR=../Rcpp/release
cd ../../multiStateFlowReliability
mkdir release
cd release
make .. -DCMAKE_BUILD_TYPE=Release  -DmaxFlowAlgorithms_DIR=../../maxFlowAlgorithms/release -DRcpp_DIR=../../Rcpp/releasemake && make install
make && make install
```