# ConnectOR-optimized
Optimized and modified version of ConnectOR (https://github.com/Carlospq/ConnectOR) originaly developed by Carlos Pulido and Daniel Kużnicki.

Optimized & modified by: Sasti Gopal Das
![Human JPG](https://github.com/cobRNA/ConnectOR-optimized/blob/main/image/human-png-34015.jpg)
## Prerequisites:
#### Python3: (install Python 3)

```sh
install Packages 

pip3 install wget
pip3 install pandas
```
#### R:
```sh
install Packages

ggplot2
tidyr
dplyr
```
#### install BedTools
```sh
sudo apt-get install bedtools
```
## Set Up:
```sh
To execute connector clone this repository locally. Following folders/files must be in the same path from command is executed:
    • ConnectOR.v.3.0.py 
    • config 
    • dictionaries 
    • ./scripts/
```
#### config example (Human vs Mouse)
```sh
species	assembly_version	annotation	chainmap
Human	hg38
Mouse	mm39
```
#### config example (Tab delimited file):
```sh
species	assembly_version	annotation	chainmap
Human	hg38	GTFs/gencode.v47.primary_assembly.annotation.gtf.gz	chainmaps/hg38ToMm39.over.chain.gz
Mouse	mm39	GTFs/gencode.vM36.primary_assembly.annotation.gtf.gz	chainmaps/mm39ToHg38.over.chain.gz
```
#### Execution command:
```sh
Python3 ConnectOR.v.3.0.py

usage: ConnectOR.v.3.0.py [-h] [-mM MINMATCH] [-g]

optional arguments:

-h, --help show this help message and exit
-mM MINMATCH, --minMatch MINMATCH
       0.N Minimum ratio of bases that must remap in liftOver
       step.Default: 30 (0.30 minimum ratio)
-g, --gene Generate results at gene level(default: False)

```
##  Gene Clustering
![Cluster JPG](https://github.com/cobRNA/ConnectOR-optimized/blob/main/image/Clustering.jpg)

##  Orthology - Human (v47) vs Mouse (vM36) - GENCODE
![othology JPG](https://github.com/cobRNA/ConnectOR-optimized/blob/main/image/othology.jpg)

