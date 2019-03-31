#!/bin/bash

./gpu -n 500 
./gpu -n 1000 
./gpu -n 2000 
./gpu -n 4000 
./gpu -n 8000 
./serial -n 500 
./serial -n 1000 
./serial -n 2000 
./serial -n 4000 
./serial -n 8000  
