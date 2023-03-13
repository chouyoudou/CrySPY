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
| Vacuum    | Float     | 0       | The vacuum layer above the upper atoms, set to avoid the effects caused by periodic conditions.           |
| thickness | Float     |         | The thickness of the upper structure, which is the thickness of the structure you want to predict by CrySPY. |


## License
CrySPY is distributed under the MIT License.  
Copyright (c) 2018 CrySPY Development Team
