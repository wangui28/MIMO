This simulation code package is mainly used to reproduce the results of the following paper [1]:

[1] C. Hu, L. Dai, T. Mir, Z. Gao, and J. Fang, "Super-resolution channel estimation for mmWave massive MIMO with hybrid precoding," IEEE Trans. Veh. Technol., vol. 67, no. 9, pp. 8954-8958, Sep. 2018.

*********************************************************************************************************************************
If you use this simulation code package in any way, please cite the original paper [1] above. 
 
The author in charge of this simulation code pacakge is: Chen Hu (email: huc16@mails.tsinghua.edu.cn).

Reference: We highly respect reproducible research, so we try to provide the simulation codes for our published papers (more information can be found at: 
http://oa.ee.tsinghua.edu.cn/dailinglong/publications/publications.html). 

Please note that the MATLAB R2012a is used for this simulation code package,  and there may be some imcompatibility problems among different MATLAB versions. 

Copyright reserved by the Broadband Communications and Signal Processing Laboratory (led by Dr. Linglong Dai), Tsinghua National Laboratory
for Information Science and Technology (TNList), Department of Electronic Engineering, Tsinghua University, Beijing 100084, China. 

*********************************************************************************************************************************
Abstract of the paper: 

Channel estimation is challenging for millimeter-wave massive
MIMO with hybrid precoding, since the number of radio frequency
chains is much smaller than that of antennas. Conventional compressive
sensing based channel estimation schemes suffer from severe resolution
loss due to the channel angle quantization. To improve the channel estimation
accuracy, we propose an iterative reweight-based superresolution
channel estimation scheme in this paper. By optimizing an objective function
through the gradient descent method, the proposed scheme can iteratively
move the estimated angle of arrivals/departures towards the optimal
solutions, and finally realize the superresolution channel estimation. In
the optimization, a weight parameter is used to control the tradeoff between
the sparsity and the data fitting error. In addition, a singular value
decomposition-based preconditioning is developed to reduce the computational
complexity of the proposed scheme. Simulation results verify the
better performance of the proposed scheme than conventional solutions.

*********************************************************************************************************************************
How to use this simulation code package?

main.m (or main_UPA.m) can figure out the NMSE curve and the average spetral efficiency curve of the proposed IR-based super-resolution channel estimation

IR_SURE_CE.m (abbr. for Iterative Reweight-based SUper-REsolution Channel Estimation) is a matlab function:
	
	OUTPUTS:
		theta_es:	estimated AoAs and AoDs
		z_es:		estimated path gains
		err:		the residue Frobenius norm

	INPUTS:
		Y:	received pilots
		X:	transmitted pilots with precoding
		W:	receiver combining matrix
		Nx,Nt,Nr,Ny:	the dimensions of the matrices
		Rth:	a residue threshold based on noise power priori. err < Rth means a successful channel estimation 
*********************************************************************************************************************************

Enjoy the reproducible research!








