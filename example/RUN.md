# Try it out

If you managed to install the library and the code you should be able to
run the test program.
The great advantage of the filter diagonalization approach
is the low storage cost when the method is coupled with a program like
the direct-ADC(4) which recalculates the matrix vector
product at each iteration.
Here, for test purposes, the symmetric randomly generated matrix is
stored in memory.

There are currently three slightly different versions of the testing program.
The `test_flset` and `test_lanczos_first` start both with a lanczos run which
generates estimation of the eigenvalues. While the `test_lanczos_first`
generates a very dense spectrum, the `test_flset` creates some spectral
regions which contain less eigenvalues. In these regions the method
should produce better results.
You will be asked for some input parameters such as the energy window and
the width of the filter function. After that you should see something like the
text below:

```
Running filter diagonalization unit test
Matrix size: 800
Enter subspace dimension L [default 8]: 23
Enter target energy window (emin emax) inside(default) -2.277302 4.322526: -0.5 1
Enter filter width sigma [suggested 0.225000]: 0.1
Using L=23, sigma=0.100000, emin=-0.500000, emax=1.000000

=========================================
| Filter Diagonalizer    C. Villani     |
|                                       |
| Last Revision   Jan 2002              |
=========================================


Memory requested for the matrix construction: 6400 bytes


------------------------------------------------------
 Matrix dimension            :   800
 Spectral boundaries         :   [-2.277 : 4.323]
 Number of basis vectors     :   23
 Width of the filter function:   0.100000
 Energy interval             :   [-0.500 : 1.000]
------------------------------------------------------

the filter diagonalizer will need at most 179096 bytes



letto input

Starting Chebishev expansion...
```

If the subspace generation has been successfull, the test program tries to
calculate the eigenvectors and estimate the residual. You should see
something like this:

```
Iteration n. 67 completed
Iteration n. 68 completed

Writing Hpsi vectors on blkfil_hpsimat file...

Writing psi transposed vectors on fil_psi file...

Reading psi transposed vectors from fil_psi file...

Filter diagonalization completed.
----------------------------
Calling analyze...

set_robin_input: ndim=800, L=23
Effective subspace dimension: 23 / 23

--- Ritz residuals ---
E[ 0] = -19.532490047045   residual = 9.056e+01
E[ 1] = -0.990187754082   residual = 8.491e+01
E[ 2] = -0.766333737768   residual = 8.419e+01
E[ 3] = -0.593195813348   residual = 8.173e+01
E[ 4] = -0.469176021390   residual = 8.757e+01
E[ 5] = -0.367060562665   residual = 8.432e+01
E[ 6] = -0.329052653209   residual = 8.122e+01
E[ 7] = -0.310135321562   residual = 8.275e+01
E[ 8] = -0.209092836552   residual = 8.345e+01
E[ 9] = -0.135389557740   residual = 8.220e+01
E[10] = -0.059389045039   residual = 8.451e+01
E[11] =  0.016515145085   residual = 8.271e+01
E[12] =  0.048802179603   residual = 8.063e+01
E[13] =  0.157558271964   residual = 8.415e+01
E[14] =  0.194300209065   residual = 8.775e+01
E[15] =  0.243113750393   residual = 7.901e+01
E[16] =  0.371327681484   residual = 8.142e+01
E[17] =  0.443027770665   residual = 8.185e+01
E[18] =  0.560465927054   residual = 8.326e+01
E[19] =  0.568961679373   residual = 8.277e+01
E[20] =  0.658116824069   residual = 8.251e+01
E[21] =  0.957762522575   residual = 8.508e+01
E[22] =  14.264364433619   residual = 9.164e+01

--- Physical eigenvalues ---
evalH[0] = -19.532490047045
evalH[1] = -0.990187754082
evalH[2] = -0.766333737768
evalH[3] = -0.593195813348
evalH[4] = -0.469176021390
evalH[5] = -0.367060562665
evalH[6] = -0.329052653209
evalH[7] = -0.310135321562
evalH[8] = -0.209092836552
evalH[9] = -0.135389557740
evalH[10] = -0.059389045039
evalH[11] = 0.016515145085
evalH[12] = 0.048802179603
evalH[13] = 0.157558271964
evalH[14] = 0.194300209065
evalH[15] = 0.243113750393
evalH[16] = 0.371327681484
evalH[17] = 0.443027770665
evalH[18] = 0.560465927054
evalH[19] = 0.568961679373
evalH[20] = 0.658116824069
evalH[21] = 0.957762522575
evalH[22] = 14.264364433619

test finished with exit code:0
```
