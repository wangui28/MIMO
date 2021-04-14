This simulation code package is mainly used to reproduce the results of the following paper [1]:

[1] X. Wei, C. Hu, and L. Dai, "Deep learning for beamspace channel estimation in millimeter-wave massive MIMO systems," IEEE Trans. Commun., vol. 69, no. 1, pp. 182-193, Jan. 2021.

*********************************************************************************************************************************
If you use this simulation code package in any way, please cite the original paper [1] above. 
 
The author in charge of this simulation code pacakge is: Xiuhong Wei (email: weixh19@mails.tsinghua.edu.cn).

Ackonwledge: We are very grateful for the author of following reference paper. Our source code is improved based on their source code. 

[29] M. Borgerding, P. Schniter, and S. Rangan, “AMP-inspired deep networks for sparse linear inverse problems,” IEEE Trans. Signal Process., vol. 65, no. 16, pp. 4293–4308, Aug. 2017.

Reference: We highly respect reproducible research, so we try to provide the simulation codes for our published papers (more information can be found at: 
http://oa.ee.tsinghua.edu.cn/dailinglong/publications/publications.html). 

Please note that the MATLAB R2020a, Python 3.6 and Tensorflow 1.10.0 are used for this simulation code package,  and there may be some imcompatibility problems among different versions. 

Copyright reserved by the Broadband Communications and Signal Processing Laboratory (led by Dr. Linglong Dai), Beijing National Research Center for Information Science and Technology (BNRist), Department of Electronic Engineering, Tsinghua University, Beijing 100084, China. 

*********************************************************************************************************************************
Abstract of the paper: 

Millimeter-wave massive multiple-input multiple-output (MIMO) can use a lens antenna array to considerably reduce the number of radio frequency (RF) chains, but channel estimation is challenging due to the number of RF chains is much smaller than that of antennas. By exploiting the sparsity of beamspace channels, the beamspace channel estimation can be formulated as a sparse signal recovery problem, which can be solved by the classical iterative algorithm named approximate message passing (AMP), and its corresponding version learned AMP (LAMP) realized by a deep neural network (DNN). However, these existing schemes cannot achieve satisfactory estimation accuracy. To improve the channel estimation performance, we propose a prior-aided Gaussian mixture LAMP (GM-LAMP) based beamspace channel estimation scheme. Specifically, based on the prior information that beamspace channel elements can be modeled by the Gaussian mixture distribution, we first derive a new shrinkage function to refine the AMP algorithm. Then, by replacing the original shrinkage function in the LAMP network with the derived Gaussian mixture shrinkage function, a prior-aided GM-LAMP network is developed to estimate the beamspace channel. Simulation results on both the theoretical channel model and the ray-tracing based channel dataset show that the proposed GM-LAMP network can achieve better estimation accuracy.


*********************************************************************************************************************************
How to use this simulation code package?

This package consists of four folders: (1) CE for SV Channel Model; (2) CE for DeepMIMO Dataset; (3) Other Results; (4) Train Network.  

(1) Open folder "CE for SV Channel Model" and run "Main.m", you will see Fig. 4 NMSE performance comparison for ULAs based on the Saleh-Valenzuela channel model (by setting type=1 in codes) and Fig. 5 NMSE performance comparison for UPAs based on the Saleh-Valenzuela channel model (by setting type=2 in codes) in our paper. 

(2) Open folder "CE for DeepMIMO" and run "Main.m", you will see Fig. 6 NMSE performance comparison for ULAs based on the DeepMIMO dataset (by setting type=1 in codes) and Fig. 7 NMSE performance comparison for UPAs based on the DeepMIMO dataset (by setting type=2 in codes) in our paper. 

(3) Open folder "Other Results" and run "Main_Select_SNR_Ranges.m", you will see Fig. 3 NMSE performance comparison of different trained GM-LAMP networks with different training settings in our paper.

(4) Open folder "Other Results" and run "Main_Complexity.m", you will see Fig. 8 the number of complex multiplications against the number of antennas N in our paper.

(5)  Open folder "Other Results" and run "Main_Sumrate.m", you will see Fig. 9 sum-rate for beam selection against different NMSE for the beamspace channel estimation in our paper.

(6) Open folder "Other Results" , then open folder "Convergence",  and run "Main_Convergence.m", you will see Fig. 10 NMSE performance against the number of layers for the GM-LAMP network in our paper.

The above results are based on our  trained networks. If you want to train new networks, please  refer to folder "Train Network".  In order to train new networks, you have to generate beamspace channel data 'x', the sensing matrix 'A', and the corresponding measurement data 'y'.  

(7) Run "generate_channel_data_for_SV_channel model.m", you can generate the beamspace channel dataset based on the Saleh-Valenzuela channel model. It is noted that the channel dataset from DeepMIMO can be obtained by the following link: DeepMIMO Dataset. [Online]. Available: http://www.DeepMIMO.net. The beamspce channel dataset corresponding with DeepMIMO can be further obtained by referring to "generate_channel_data_for_SV_channel model.m".

(8) Run "generate_CS.m", you can generate a new compressive sensing matrix. It is noted that different compressive sensing matrices need to train different networks.

(9) Run "generate_measurement.py" by setting the correct file name for loading channel data and compressive sensing matrix according to the generated files in step (7) and step (8) and setting the expected SNR ranges, you can generate the corresponding measurement dataset.  

(10) Run "Train_GM-LAMP.py"  by setting the correct file name for loading channel data, compressive sensing matrix and measurement data according to the generated files in step (7), step (8), and step (9),  you will obtain a .npz file for saving the trained network. 

(11) Run "load_trained_network.py" by setting the corret .npz file name for loading the trained network according to generated files in step (10), you will obtain a .mat file for saving the network's parameters. The generated .mat file can be used to reconfigure the parameters in "GM-LAMP.m".

It is noted that there may be some differences in the results of different training processes. 

*********************************************************************************************************************************
Enjoy the reproducible research!