##Models tested in ModelAssess
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

Modelassess will look for configuration files $HOME/.modelassess.ini and
$PWD/modelassess.ini.  Configurations in the second file will be laid over the
first (i.e. models turned off in $HOME/.modelassess.ini can be switched back on
in $PWD/modelassess.ini). The format is standard INI, and models are specified
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
