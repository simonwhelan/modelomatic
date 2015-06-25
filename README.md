#Modelomatic
---
_Written by Simon Whelan, with contributions from James Allen, Ben Blackburne and David Talavera_



ModelOMatic is a C++ program designed for rapid phylogenetic model selection on protein coding genes. Please see the [manual](https://drive.google.com/file/d/0B2HVW-VOuzH-MjRKaVcwd0VPN1E/edit?usp=sharing) for details of the program and its settings.

ModelOMatic is distributed free of charge for personal and academic use only under a [GNU GPL v3](http://www.gnu.org/licenses/gpl.html). Modelomatic is a piece of academically developed software and should be used at your own risk! 

Pay attention to the output and use your own judgements when performing analysis. Modelomatic uses lots of heuristics that seem to perform well when tested across the PANDIT database, but they may not perform well on your data. 

##Download information
* Precompiled binaries are available [here](https://drive.google.com/folderview?id=0B2HVW-VOuzH-U0NsV0hWMkxEQTQ&usp=sharing)
* The user manual is available [here](https://drive.google.com/file/d/0B2HVW-VOuzH-MjRKaVcwd0VPN1E/edit?usp=sharing)
* Example files are available [here](https://drive.google.com/folderview?id=0B2HVW-VOuzH-SlVaRzM1bkl0cDA&usp=sharing)

##Version history
Modelomatic will be periodically updated and changes that may affect the outcome of your analyses will be announced here. Minor fixes that address bugs that are not expected to change results are not included here.  

##Models tested in Modelomatic
* RY
* DNA
    - JC
    - FEL
    - K2P
    - HKY
    - REV 
* AA 
    - EQU
    - BLOSUM62
    - DAY
    - JTT
    - WAG 
    - LG
    - VT
    - cpREV
    - mtArt
    - mtMam
    - mtREV
    - rtREV
    - HIVb
    - HIVw 
* Codon
    - F0
    - F1X4
    - F3X4
    - F64

##How to specify models with a configuration file

Modelomatic will look for configuration files $HOME/.modelomatic.ini and
$PWD/modelomatic.ini.  Configurations in the second file will be laid over the
first (i.e. models turned off in $HOME/.modelomatic.ini can be switched back on
in $PWD/modelomatic.ini). The format is standard INI, and models are specified
with 0/1, as below:

~~~~~~

[Amino Acid]
;Ignore Viral models
rtREV=0
HIVb=0
HIVw=0

[Codon]
F0=0

~~~~~~

Currently only Amino Acid and Codon models can be specified. The default case
is to run everything. The INI file is case insensitive.
