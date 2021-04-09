# Computing Torsion Angles for N-glycan

Understanding common torsion angles in N-glycans is important in the correct modelling of structures from experimental data. These scripts allow torsion angles data collection to be automated using [Privateer](http://legacy.ccp4.ac.uk/html/privateer.html).

![Torsion Angles](https://github.com/Dialpuri/N-glycanTorsionAngles/blob/269acd1259f91895f2b1bf4646da6b126afe55d8/TorsionAngles.png)
(Image From: https://doi.org/10.1107/s2059798316016910)
## Data
The data, both raw and calculated is available from [here](https://drive.google.com/drive/folders/1vUMD7k50AcrIeOBhN2oZg_b-LH9Ea-qB?usp=sharing)

## Version History

### Possible Sugar Connections
* Took in an input file of computed data from privateer and looked for consecutive positions in the carbohydrate chain
* Printed this output for manual visualisation

### Auto Version 1.0
* Extension of Possible Sugar Connections, consecutive sugars identified
* Downloaded the PDB file using Biopython
* Computed the torsion angles between carbohydrates and between residue and carbohydrate using Privateer's print_glycosidic_torsions function
* Outputted file with torsion angle data in json format

### Auto Version 1.1
* Extension of auto_v1.0 
* Added correct RSCC for both sugars involved 
* Residue-carbohydrate output fixed to only be for RSCC > 0.8

### Auto Version 1.2
* Using custom function in Privateer.py to gather all torsion angles for every glycan
* Outputted only those which were in the input file (i.e. RSCC > 0.8)
* Run time reduced significantly

### Auto Version 1.3
* Fixed error in extra data in output file
* Implemented custom function to get residue-carbohydrate from Privateer
* Run time reduced significantly

### Auto Version 1.3.1
* Fixed error in output header

## Installation

To run the latest script (1.3.1), you will need the following libraries installed:
* Privateer
  * Installation instructions [here](https://github.com/glycojones/privateer)
  * Or included as part of [CCP4](https://www.ccp4.ac.uk/)
* Biopython (Version <= 1.7.6)
* Pandas

using the following code in the CCP4 console:
```
pip install biopython
pip install pandas
```
Replace the privateer.py file in C:\CCP4-7\7.1\Lib\py2 with the one here

## Usage

From the CCP4-Console run:
```
ccp4-python auto_v1.3.1.py
```

## Explanation 
See Code Explanation.pdf for informal explanation of the code

## Data Splitting
The initial data was organised by RSCC > 0.8 and diagnostically yes and split into 20 chunks using files inside the Data Splitting folder
