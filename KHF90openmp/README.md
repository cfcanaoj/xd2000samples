# Kelvin-Helmholtz Instability Test
This code solves hydrodynamics of  Kelvin-Helmholtz Instability.

## Description of the problem

https://www.astro.princeton.edu/~jstone/Athena/tests/kh/kh.html

## Numerical setup

https://ui.adsabs.harvard.edu/abs/2012ApJS..201...18M/abstract

# How to run

## How to copy the source code
After you login the server, `xd01.cfca.nao.ac.jp`, follow the instruction.

    cd /work/<username>
    git clone git@github.com:cfcanaoj/xd2000samples test
    
## Change directory
You can change the directory by `cd`. `test` can be different if you change the name above.
    
    cd test/KHF90openmp

## Compile code
First you need to prepare Intel one-api environment. This setup is also need to run the program. You should write it in `.basrhrc`.
    
    source enable-oneapi.sh
    
After that you can compile it. If you want to change the compile option, edit `Makefile`.
    
    make
    
 ## Run program
 Before the submission, you need to edit the batch script.
    
    vim slm_xd.sh

In the following part, you need to specify your own partition.

    #SBATCH --partition=M-test-cfca

Then you can submit the job.

    sbatch slm_xc.sh
    
You can confirm the job by the following command.
    
    squeue --me
   

# How to confirm the results
Let us move to analysis server.

    ssh an10@cfca.nao.ac.jp
    cp /xd-work/<username>/test/KHF90openmp .
    cd KHF90openmp
    module load gnuplot
    python MakePlot.py
    display images/den00100.png


