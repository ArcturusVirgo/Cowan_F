cmake_minimum_required(VERSION 3.25)
project(Cowan Fortran)
enable_language(Fortran)


add_executable(widen_4 widen_4.f90 module/DFile_Mod.f90 module/cowan_data_processing.f90)
add_executable(fraction fraction.f90 module/DFile_Mod.f90)

add_executable(test test.f90 module/DFile_Mod.f90 module/ion_abundance.f90)
add_executable(SIM SIM.f90)