# ConnectOR-optimized
Optimized and modified version of ConnectOR (https://github.com/Carlospq/ConnectOR) originaly developed by Carlos Pulido and Daniel Kużnicki.

Optimized & modified by: Sasti Gopal Das
![Human JPG](https://github.com/cobRNA/ConnectOR-optimized/blob/main/raw/ConnectOR.jpg)
## Prerequisites:
#### Python3:

```sh
install Packages

wget
pandas
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
    • ConnectOR_V2.1.py 
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
#### config example (Human vs Zebrafish)
```sh
species	assembly_version	annotation	chainmap
Human	hg38
Zebrafish	danrer11
```
#### Execution command:
```sh
Python ConnectOR_v2.1.py

usage: ConnectOR.v2.py [-h] [-mM MINMATCH] [-g]

optional arguments:

-h, --help show this help message and exit
-mM MINMATCH, --minMatch MINMATCH
       0.N Minimum ratio of bases that must remap in liftOver
       step.Default: 30 (0.30 minimum ratio)
-g, --gene Generate results at gene level results (default: False)

```
##  Gene Clustering
![Cluster JPG](https://github.com/cobRNA/ConnectOR-optimized/blob/main/raw/Cluster.jpg)

