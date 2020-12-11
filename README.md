# ElasticGen

Basic Elastic generator using [TGenPhaseSpace](https://root.cern.ch/doc/master/classTGenPhaseSpace.html) to generate events and output to LUND format for input into the clas12 simulation chain. 

##### Build ElasticGen:

```shell
git clone https://github.com/uofscphysics/ElasticGen.git
mkdir -p ElasticGen/build
cd ElasticGen/build
cmake .. 
make
```

##### Run ElasticGen:

```shell
# ElasticGen | outputfile | num of event | beam energy | min Q2 | max Q2
./ElasticGen outputfile.lund 10000 10.6 0.0 20
```

##### TODO [If I need to or if there is time]:

- [ ] Add Radiative Corrections option

- [ ] Add options as command line arguments with [clipp](https://github.com/muellan/clipp)

- [ ] Add cross section information into generator
