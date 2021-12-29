# Fast explicit dynamics finite element algorithm for transient heat transfer (MIT License)

![fig1](https://user-images.githubusercontent.com/93865598/147568282-8c4247c0-cab2-4636-8dda-227df9d0f58c.PNG)
This is the source repository for the paper:

Zhang, J., & Chauhan, S. (2019). Fast explicit dynamics finite element algorithm for transient heat transfer. International Journal of Thermal Sciences, 139, 160-175. [doi:10.1016/j.ijthermalsci.2019.01.030](https://www.sciencedirect.com/science/article/abs/pii/S1290072918317186).

Please cite the above paper if you use this code for your research.

# Environment:
•	Windows 10

•	Visual Studio 2017

•	OpenMP
# How to build:
1.	Download the source repository.
2.	Visual Studio 2017->Create New Project (Empty Project)->Project->Add Existing Item->FEDFEM.cpp.
3.	Project->Properties->C/C++->Language->OpenMP Support->Yes (/openmp).
4.	Build Solution (Release/x64).
# How to use:
1.	(cmd)Command Prompt-> ![fig2](https://user-images.githubusercontent.com/93865598/147568308-8752fdbb-8067-4c3d-9089-13f631476ce4.PNG)
2.	Output: T.vtk
# How to visualize:
1.	Open T.vtk. (such as using ParaView)

![fig3](https://user-images.githubusercontent.com/93865598/147568315-1d2c3f4c-4dd7-4a6e-b3fe-4169c26555c7.PNG)
# Boundary condition (BC):
1.	Node index: HFlux, Convc, Radia, FixT.
2.	Element index: BodyHFlux.
# Notes:
1.	Node and Element index can start at 0, 1, or any but must be consistent in a file.
2.	T_Iso.txt, T_Ortho.txt, T_Aniso.txt, T_Convc.txt, T_Radia.txt, T_BodyHFlux.txt, node and element index start at 0.
3.	T_Iso_n1.txt, T_Ortho_n1.txt, T_Aniso_n1.txt, T_Convc_n1.txt, T_Radia_n1.txt, T_BodyHFlux_n1.txt, node and element index start at 1.
