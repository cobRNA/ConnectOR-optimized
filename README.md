# ConnectOR-optimized
Optimized and modified version of ConnectOR (https://github.com/Carlospq/ConnectOR) originaly developed by Carlos Pulido and Daniel Kużnicki.

Optimized & modified by: Sasti Gopal Das
## Prerequisites:
#### Python3:
```sh
python -m pip install wget
python -m pip install pandas
```
### install BedTools
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
```




