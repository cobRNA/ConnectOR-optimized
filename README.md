# ConnectOR-optimized
Optimized and modified version of ConnectOR (https://github.com/Carlospq/ConnectOR) originaly developed by Carlos Pulido and Daniel Kużnicki.

Optimized & modified by: Sasti Gopal Das
## Prerequisites:
#### install python library
```sh
pip install argparse
pip install wget
pip install pandas.
pip install matplotlib.
pip install networkx.
pip install matplotlib-venn
```
#### install python library
```sh
install.packages("reshape2")
install.packages("ggplot2")
install.packages("scales")
install.packages("dplyr")
```
### install BedTools
```sh
BedTools
or
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
#### Orthology - Human (v43) vs Mouse (vM32) – GENCODE

### Notes:
1. **Ensure the Image is in the Correct Path**: The `Human_Mouse.png` file should be in the root directory of your repository, or adjust the path accordingly if it's in a subdirectory.
2. **Correct Markdown Syntax**: The correct syntax to display an image in Markdown is `![Description of image](path/to/image.png)`. Ensure no extra brackets or parentheses are used.
3. **Check Image Visibility**: After pushing the changes, verify on GitHub that the image is displayed correctly in the `README.md` file.

### Full Sequence of Commands

Here’s the full sequence of commands for uploading the image and updating the `README.md` file:

```sh
# Navigate to your local repository
cd (https://github.com/cobRNA/ConnectOR-optimized/edit/main/)



# Add and commit the image file
git add Human_Mouse.png
git commit -m "Add Human_Mouse.png"
git push origin main

# Update the README.md file to include the image
echo '![Description of image](Human_Mouse.png)' >> README.md

# Add and commit the updated README.md file
git add README.md
git commit -m "Update README to include image"
git push origin main

