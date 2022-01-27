# Fast explicit dynamics finite element algorithm for transient heat transfer (MIT License)
![GitHub](https://img.shields.io/github/license/jinaojakezhang/fedfem)
![GitHub top language](https://img.shields.io/github/languages/top/jinaojakezhang/FEDFEM)
<p align="center"><img src="https://user-images.githubusercontent.com/93865598/147568282-8c4247c0-cab2-4636-8dda-227df9d0f58c.PNG"></p>
This is the source repository for the paper:

| Zhang, J., & Chauhan, S. (2019). Fast explicit dynamics finite element algorithm for transient heat transfer. *International Journal of Thermal Sciences*, 139, 160-175. [doi:10.1016/j.ijthermalsci.2019.01.030](https://www.sciencedirect.com/science/article/abs/pii/S1290072918317186). |
| --- |

Please cite the above paper if you use this code for your research.

If this code is helpful in your projects, please help to :star: this repo or recommend it to your friends. Thanks:blush:
## Environment:
- Windows 10
- Visual Studio 2017
-	OpenMP
## How to build:
1.	Download the source repository.
2.	Visual Studio 2017->Create New Project (Empty Project)->Project->Add Existing Item->FEDFEM.cpp.
3.	Project->Properties->C/C++->Language->OpenMP Support->**Yes (/openmp)**.
4.	Build Solution (Release/x64).
## How to use:
1.	(cmd)Command Prompt->build path>project_name.exe input.txt. Example: <p align="center"><img src="https://user-images.githubusercontent.com/93865598/147874639-c1cecfa9-8405-4006-b907-4e5751491dd4.PNG"></p>
2.	Output: T.vtk
## How to visualize:
1.	Open T.vtk. (such as using ParaView)
<p align="center"><img src="https://user-images.githubusercontent.com/93865598/147568315-1d2c3f4c-4dd7-4a6e-b3fe-4169c26555c7.PNG"></p>

## How to make input.txt:
1.	T_Iso.inp (Abaqus input) is provided in the “models”, which was used to create T_Iso_n1.txt.
## Material types:
1.	Isotropic, orthotropic, and anisotropic thermal conductivities.
## Boundary conditions (BCs):
1.	Node index: HFlux, Convc, Radia, FixT.
2.	Element index: BodyHFlux.
## Notes:
1.	Node and Element index can start at 0, 1, or any but must be consistent in a file.
2.	Index starts at 0: *.txt.
3.	Index starts at 1: *_n1.txt.
## Feedback:
Please send an email to jinao.zhang@hotmail.com. Thanks for your valuable feedback and suggestions.
