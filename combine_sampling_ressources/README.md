# Combine samplings (shots) from multiple runs of on quantum computers
This directory's contents can be copied in your dedicated /run folder on your computer. If you follow the following folder format for your runs, this code will take care of combining the matrices for you.

## Format
Here is how you need to name your directories containing your runs:`{N}sites_mu{mu}_{device}_date_month_{#_of_run}_{shots}sh/`.

For example, you could run:
```bash
mkdir 2sites_mu2_sherbrooke_27_juil_2_8000sh/
```

Then, in this file, you would have your parameters.py file, as well as your excitation.def file, if applicable. 

## Execution
To run the code:
```bash
python3 sampling.py
```

and follow the instructions.
