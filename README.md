# qPCR_analysis

The sets of code provided here are to anaylize data generated from RT-QPCR.
This will take a CSV file with sample names and gene names and calculate the fold changes in mRNA expression. 
(Tested on linux only, but should work on MACOS)

To install run:
git clone https://github.com/sudo-shaka/qPCR_analysis.git
cd qPCR_anaylisis && make
To make executable from anywhere for to path: example: /usr/bin /usr/local/bin or .local/bin
type ./qPCR /path/to/csv/to/anaylize/ and it will calculate your fold change

Alternatively there is an python script for those more famliar with that language
