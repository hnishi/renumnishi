# renumnishi
re-numbering all residue in pdb file

## BUILD

edit "Makefile" and then

```
$ make clean
$ make
```

## USAGE

create input file for program.   
see inp_template.txt   
you have to fill three parameters;   
i.e.   
INPUTPDB1 (input pdb file to be renumbered)   
OUTPUTFILE1 (renumbered pdb file (output))   
START_RES (starting residue number (interger))   
```
$ ./a.out inp.txt
```

