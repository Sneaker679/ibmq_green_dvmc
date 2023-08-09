# Combine samplings (shots) from multiple runs of on quantum computers
This directory's contents can be copied in your dedicated /run folder on your computer. If you follow the following folder format for your runs, this code will (hopefully) take care of combining the matrices for you.

## Format
Here is how you need to name your directories containing your runs:`{N}sites_mu{mu}_{device}_{date}_{month}_{#_of_run}_{shots}sh/`.

For example, you could run:
```bash
mkdir 2sites_mu2_sherbrooke_27_juil_2_8000sh/
```

Then, in this file, you would have your parameters.py file, as well as your excitation.def file, if applicable. You then run the simulation like usual. Then, to use the combine the results with this code: 

## Execution
To run the code:
```bash
python3 sampling.py
```
if you want to combine the results with an average. Or you can:
```bash
python3 median.py
```
if you want to combine the results with a median.

Simply follow the instructions on the screen.
