![cryspy_logo](./cryspy_fix-03.png)

# CrySPY
CrySPY is a crystal structure prediction tool written in Python.  
Document site is moved to https://tomoki-yamashita.github.io/CrySPY_doc/

## Important changes
- Now CrySPY can be used for interface structure prediction.
- New interface for Atomic Simulation Environment (ASE).


## Parameters of interface structure prediction
| Name      | Value     | Default | Description                                                                                               |
| --------- | --------- | ------- | --------------------------------------------------------------------------------------------------------- |
| buffer    | Float     | 0       | Distance between substrate and upper layer, set to prevent overlap of the upper atoms and substrate atoms of the initial structure. |
| vacuum    | Float     | 0       | The vacuum layer above the upper atoms, set to avoid the effects caused by periodic conditions.           |
| thickness | Float     |         | The thickness of the upper structure, which is the thickness of the structure you want to predict by CrySPY. |
| up_atype | Str     |         | Atomic element of the upper structure. e.g. up_atype = Al O |
| up_nat | int     |         | Number of atoms in upper structure. |
| struc_mode | interface, crystal, mol, mol_bs     |         | Structure generation mode. |
<img width="252" alt="截屏2023-03-13 16 12 22" src="https://user-images.githubusercontent.com/60209970/224632007-b6480cd6-b27b-452c-bb5d-b94968e6a871.png">
![屏幕录制2022-08-08 15 28 40 (1)](https://user-images.githubusercontent.com/60209970/224637116-156c1b6c-79e8-4db5-8e8a-9d9bfd3141f9.gif)


## License
CrySPY is distributed under the MIT License.  
Copyright (c) 2018 CrySPY Development Team
