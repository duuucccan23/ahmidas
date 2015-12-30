# Introduction #

Because of the unforeseeable diversity of systems that ahmidas may run on, it is important to begin devising a set of unit tests as early on in the development process as possible. Many core components should be rigorously tested and testable.

# Tests #

## Endianness ##

A program that:
  1. should test and report on the endianness of the system ahmidas is running on (or going to run on in case of cross-compilation)
  1. verifies proper conversion of little endian to big endian and vice versa
  1. does the above for at least double and float, but perhaps more

## SU3 ##

### Matrix ###

#### Reunitarization ####

A program that:
  1. creates a random (full rank) matrix and then reunitarizes it and verifies that the result is indeed special unitary
  1. creates a random SU3 matrix and verifies that it is unchanged under reunitarization

### Vector ###

## MPI ##

## Grid ##

## Fields ##

### Gauge ###

A program that:
  1. Reads a complete gauge field (should be of reasonable size like 4^4)
  1. Writes a complete gauge field

Perhaps it is useful to also test the entire reading, endianness swapping and float to double conversion/reunitarization chain here as well

#### APE smearing ####

### Spinor ###