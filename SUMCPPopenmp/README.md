# Brief description
This code just take a sum of 2D array.
```math
 s = \sum_{i,j} S_{i,j}
```
# Download
Download this at the login node of CfCA server.

	ssh <username>@xd01.cfca.nao.ac.jp
	cd /work/<username>
	git clone git@github.com:cfcanaoj/xd2000samples test
	cd test/SUMCPPopenmp

# Compile code
First you need to prepare Intel one-api environment. This setup is also need to run the program. You should write it in `.basrhrc`.
    
    source enable-oneapi.sh
    
 After that you can compile it. If you want to change the compile option, edit `Makefile`.
    
    make

 # Run program
 Before the submission, you need to edit the batch script.
    
    vim slm_xd.sh

In the following part, you need to specify your own partition.

    #SBATCH --partition=M-test-cfca

Then you can submit the job.

    sbatch slm_xc.sh
    
You can confirm the job by the following command.
    
    squeue --me
   
