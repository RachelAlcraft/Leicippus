﻿# CMakeList.txt : CMake project for Leucippus, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)

project ("Leucippus")

# Add source to this project's executable.
add_executable (Leucippus "Leucippus.cpp" "Leucippus.h" "Lattice.h" "Lattice.cpp"   "Xray.h" "Xray.cpp" "VectorThree.h" "VectorThree.cpp" "Atom.cpp" "Atom.h" "StructureFactors.h" "StructureFactors.cpp" "PeriodicTable.h" "PeriodicTable.cpp" "CifFile.h" "CifFile.cpp" "MolRep.h" "MolRep.cpp"  "Protein.h" "Protein.cpp" "ElectronDensity.h" "ElectronDensity.cpp" "Ccp4File.cpp" "Ccp4File.h" "Helper.h" "Helper.cpp" "MatrixThreeThree.h" "MatrixThreeThree.cpp" "ED_SF.cpp" "ED_SF.h")

# TODO: Add tests and install targets if needed.
