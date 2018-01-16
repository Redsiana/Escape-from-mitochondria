# Escape-from-mitochondria
The program was developped in R version 3.4.1 (2017-06-30) running on a x86_64-w64-mingw32 platform.

How to run the program to replicate or expand the results presented in the manuscript:
1. Install all required libraries by running the file libraries.R
2. Open the file PARAM_TABLE.R and choose the parameter space to explore
3. Run the script MAIN.R (warning: with a parameter space like the one used in the paper, this starts a sequence of many, heavy simulations, that can last for 2 days on one core)
4. Generate figures by running the desired Fig.x_script.R files
