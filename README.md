### Rclonte-M
**Reference: Toward An Advanced Method for Full-Waveform Hyperspectral LiDAR Data Processing**

This file is generated on 17th Jan., 2024 by its author Jie Bai (白杰). Codebase is freely available to public now.

- **Sharing/access Information——————————————————————————————**

Copyright (c) 2024 Jie Bai

Permission is hereby granted, free of charge, to any person obtaining a copy of this code and to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, and utilize it in commercial softwares.

However, note that users' improved or modified versions of this code are welcome and permitted only when acknowledging our article and with explanations in appropriate places in their articles. These improvements include, but are not limited to, the use of this code for the processing of single-wavelength waveforms and modifications to the kernel transform for advanced multispectral or hyperspectral waveform applications.

- **About The Repository——————————————————————————————————**

**Code Files:**

 measuredDataDecompositionRclonte_M.py (main part)

 RClonte_M_Method.py (method part)

**ExampleFiles**

The author stores the hyperspectral waveform samples for users to implement our code, ranging from 409 nm ~ 914 nm in this folder. Notably, the middle part of wavelengths about 491 nm ~ 800 nm owns high Signal-to-Noise Ratios (SNR) and the other redundant wavelengths have been filtered at the quality check stage in preprocessing. The users can download this hyperspectral .csv files and then change the corresponding rootpath variable in **measuredDataDecompositionRclonte_M.py (main part)** on the basis of their needs. The code will output the related decomposed results for visualization.

- **Workflow of Rclonte-M——————————————————————————————————**

![Fig  5算法流程图](https://github.com/Jie-Bai/Rclonte-M-TGRS/assets/37448920/6de861ba-9bed-42d8-8512-ddab9686f9cc)


Fig. 1. Workflow of the Rclonte-M decomposition method.






