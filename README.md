Programs for our cross-layer adaptive reference approach for SHVC, in SNR scalability at inter-mode.

# Cross-Layer Adaptive Reference based SHM Encoder
This encoder is used to evaluate the performance of our CLAR methodology [1] The main part is modified from the standard reference software SHM 12.4 [2], 
coded with C++. To work out bandwidth usage and error propagation problems resulting from information loss in time-varing networks, reference pictures 
with high quality for base layer (BL) are derived according to conditions of pictures loss at enhancement layer (EL) at deocoder end. Intuitively, the proposed 
CLAR method can save more bitrate compared with the anchor for equal quality. Moreover, this scheme eases the effect of lost pictures at decoder end.

# References
[1] Liu X, Yuan J, Xie R, et al. Enhanced inter-layer prediction for SHVC based on cross-layer adaptive reference[C]. International Forum on Digital TV and Wireless Multimedia Communications. Singapore: Springer Singapore, 2021: 412-425.
[2] JCT-VC, “SHM Software,” [Online]. Available: https://hevc.hhi.fraunhofer.de/svn/svn_SHVCSoftware/tags/SHM-12.4/.
