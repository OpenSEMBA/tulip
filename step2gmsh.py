#!/usr/bin/env python

import sys
import argparse
from src.mesher import Mesher
    
def launcher(fn):
    mesher = Mesher()
    mesher.runFromInput(fn)
    

if __name__ == '__main__':
    print("-- Launching step2gmsh")

    parser = argparse.ArgumentParser(
        prog="step2gmsh",
        description=
        "Converts step CAD file into a MFEM compatible mesh.\n"
        "Please look at README.md and LICENSE for more info at:\n"
        " https://github.com/OpenSEMBA/step2gmsh"
    )
    parser.add_argument(
        "-i",
        "--input",
        help="step input file",
        type=argparse.FileType('r'),
        required=True
    )

    args = parser.parse_args()
    launcher(args.input.name)
