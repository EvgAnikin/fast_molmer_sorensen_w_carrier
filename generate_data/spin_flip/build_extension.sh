#!/bin/bash
f2py -c spin_flip.f95 --fcompiler=gnu95 --f90flags=-ffixed-line-length-none -m spin_flip
