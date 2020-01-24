# Advance Visualisation

Visualisation of the BitEpi output in Cytoscape using a Python API.  

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

### Prerequisites

Python3.6 or above

Install Pandas

```
pip install pandas
```

Install py2cytoscape

* [py2cytoscape](https://py2cytoscape.readthedocs.io/en/latest/#installation) - A collection of utilities that enables one to use Cytoscape using Python

Install Cytoscape

* [Cytoscape](https://cytoscape.org/download.html) - Platform used for visualisation

### Installing

A step by step series of instructions that tell you how to get the development env running

Clone this repository upon installation of the above prerequisites. Open Cytoscape and have it running in the background. Then open the command prompt/ terminal and navigate to the folder where the program is by typing the following command.

```
cd BitEpi/AdvanceVisualisation/
```

Then in order to run the program make sure to give a .csv file as the input to the command line argument -i as shown below. 

```
python3 ./main.py -i out1.Alpha.4.csv
```

Then Python will automatically provide you with the Cytoscape visualisation of the content in the file provided. 

## Authors

* **Milindi Kodikara**

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## Acknowledgments

* Arash Bayat

