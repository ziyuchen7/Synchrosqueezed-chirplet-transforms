# MATLAB code for Synchrosqueezed Chirplet Transforms
MATLAB code for the implementation of the Synchrosqueezed chirplet transforms with some examples.

Reference: 
@article{chen2023disentangling,
  title={Disentangling modes with crossover instantaneous frequencies by synchrosqueezed chirplet transforms, from theory to application},
  author={Chen, Ziyu and Wu, Hau-Tieng},
  journal={Applied and Computational Harmonic Analysis},
  volume={62},
  pages={84--122},
  year={2023},
  publisher={Elsevier}
}

File explanation:

dwindow.m -- compute the derivative of a given window

sqSTFT.m -- STFT-based Synchrosqueezed transforms (SST)

sqSTCT.m -- STFT-based Synchrosqueezed chirplet transforms (SCT, including chirplet transforms)

sqSTFTbase2nd.m -- STFT-based second-order SST

SCTinverse.m -- inverse map of SCT

example1.m -- example of standard linear chirps

example2.m -- example of a signal with time-varying chirp rate

example3.m -- example of wolf howling signals

![test](./example1-CT3Dview.gif)

*3D visualization of the chirplet transform*

![test](./example1-SCT3Dview.gif)

*3D visualization of the synchrosqueezed chirplet transform*

Any questions can be sent to ziyuchen@umass.edu.
