# qf15

`qf15` is a software package accompanying paper 

> *Martin Smid. Econometrics of Zero Intelligence Models by Means of L1 Data. submitted to Quantitative Finance, 2015.* 

This file contains an information how to install the package and replicate results from the paper. 

##  Prerequisities

The following conditions must be fulfilled for `qf15` to run on Linux. Running on other OS's are analogous but may require changes in source code.


1. folders `hfd`, with all permisions present elsewhere, containinng subfolders
	* `data`
	* `latex`
	* `inter/csv` 

2. file `index.csv` of zero length present in the `hfd` folder

3. program `7z` installed somewhere in path

## Building

`qf15` is built by means of `Code::Blocks`, project `qf15.cbp`. Fro the program to be built, it must be

1. `boost` library installed somewhere on computer, its include path added to `Build Options/Search directories/Compiler' of the project

2. `nlopt` library installed somewhere, its include path added to the project (as above) and a path to the library itself added to `Build Options/Search directories/Linker' of the project

The executable is then build e.g., by `Build/Build` of `Code::Blocks`.

## Usage

### Data import

Data are imported from the format of trade and quote tada provided by [www.tickdata.com] as of year 2009, i.e., one zipfile per one day, where the file for quotes of *`yyyy.mm.dd`* are stored in *`dataroot/yyyy/mm`*`/QUOTES` where *`dataroot`* is a foldear all the data are stored. The trades find itself in *`dataroot/yyyy/mm`*`/TRADES`

The data are imported by 

> `qf15 I`  *`hfdroot dataroot`* 

where *`hfdroot`* is the folder in which `hfd` finds itself. Notice, that the file `intex.csv` has to be empty before the import. To replecate the results from the paper, data from March and April 2009 must be present in the source folder.

### Replication of the paper's results

To replicate `Appendix C1`, use

> `qf15 z` *`hfdroot`*

To replicate `Appendix C2`, use

> `qf15 g` *`hfdroot`*


## Data strucure

An alternative to importing commercial data could be to create the internal database independenty. Its structure is as follows: data of each stock-market pair are storded in a separate file in 

> *`hfdroot`*`/data/`*`stock`*`_`*`market`*`_`*`yyyymmdd`*.bin

such that the individual records are binary images if `struct hfdrecord` (for its structure, see `hfd.hpp`). 
The list of all those stockpair-days is stored in 

> *`hfdroot`*`/index.csv`

which is a csv file without heading where each line has a structure

> *`yyyy,mm,dd,stock,market,datatype`*

where *`datatype`* may be presently equal only to `1`. For rescords for each stock-market pair must be stored time ascending way.
 