### Rclonte-M
**Reference: Toward An Advanced Method for Full-Waveform Hyperspectral LiDAR Data Processing**

This file is generated by its author Jie Bai (白杰) on 17th Jan., 2024. Codebase is freely available to public now 

由于我没有在网上找到任何一段完整的波形分解的代码，所有代码均为我本人自己摸索撰写，于投稿之初主动公开，供大家参考。
代码是非常珍贵的，这是我熬了无数大夜写出来的，但代码肯定存在不足之处，还望各位专家拍砖指正（之后身体指标也出现各种异常，一度到思考放弃科研）。
现在的算法、模型文章虽然很多，但愿意公开免费给出代码的却寥寥无几。曾经礼貌咨询一位学者（非学生）是否可以学习某代码，被正常收要1000元人民币，我不舍，于是至今也没有搞明白某算法，这种情况下我也曾和我同学进行过探讨、尝试复现，但是我以为的复现在没有源代码的情况即使复现了、效果不好，在假设源代码鲁棒性很好排除数据集影响因素以外就能说明人家原作者代码本身不好并进行评论吗，私认为未必。当然，我也希望能有越来越多的学者可以公开免费源代码供拙人（单指我自己）进步！



- **Sharing/access Information——————————————————————————————**

Copyright (c) 2024 Jie Bai

Permission is hereby granted, free of charge, to any person obtaining a copy of this code in their studies by citing this paper and to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, and utilize it in commercial softwares.

Also, note that users' improved or modified versions of this code are welcome and permitted only when acknowledging our article and with explanations in appropriate places in their articles. These improvements include, but are not limited to, the use of this code for the processing of single-wavelength waveforms and modifications to the kernel transform for advanced multispectral or hyperspectral waveform applications.

- **About The Repository——————————————————————————————————**

**Code Files:**

 measuredDataDecompositionRclonte_M.py (main part)

 RClonte_M_Method.py (method part)

**ExampleFiles**

The author stores the hyperspectral waveform samples for users to implement our code, ranging from 409 nm ~ 914 nm in this folder. Notably, the middle part of wavelengths about 491 nm ~ 800 nm owns high Signal-to-Noise Ratios (SNR) and the other redundant wavelengths have been filtered at the quality check stage in preprocessing. The users can download this hyperspectral .csv files and then change the corresponding rootpath variable in **measuredDataDecompositionRclonte_M.py (main part)** on the basis of their needs. The code will output the related decomposed results for visualization.

- **Workflow of Rclonte-M——————————————————————————————————**

![Fig  5算法流程图](https://github.com/Jie-Bai/Rclonte-M-TGRS/assets/37448920/6de861ba-9bed-42d8-8512-ddab9686f9cc)


Fig. 1. Workflow of the Rclonte-M decomposition method.

### NOTE
Compared with the waveform processing algorithm that takes the median after the target center position is arranged in descending order (in this paper Bai et al., 2024, IEEE TGRS), our RSE (Bai et al., 2024, RSE) paper adopts different rule to compensate for hidden echoes or weak echoes. 

Bai J, Niu Z*, Huang Y, Wang L*, et al (2024). Full-waveform hyperspectral LiDAR data decomposition via ranking central locations of natural target echoes (Rclonte) at different wavelengths. Remote Sensing of Environment, 310:114227. link this: https://www.sciencedirect.com/science/article/abs/pii/S0034425724002451?via%3Dihub

More detailed core difference is here for the two and can help readers know more:
### 1. **Algorithm**:
The underlying idea of the algorithms is called Ranking Central Locations of Natural Target Echoes. But the steps after the descending arrangement in re-optimization are the most critical and determine the fundamental differences between the two. The Rclonte algorithm proposed by us makes detailed screening criteria after the descending order of the central position (Rclonte), including (see page 6 of the RSE article from the last 14 rows of the left column to the end of the first paragraph of the right column) 

  **a.** First, the number of all centers in the initial same order set must be greater than 1/8 times the number of wavelengths, so that the components corresponding to the initial same order can be presumed to be real at each wavelength, and most likely belong to the same component.

  **b.** However, because hidden waves and weak waves may also exist in the waveform, and the central positions of other components caused by false detection or missing detection will be carried into the current order components on this basis, so the central positions belonging to the same order may not all come from the same target sample components, we need to do further processing.

  **c.** In order to deal with the above phenomenon, the central location belonging to the same order is listed as a separate new set, and further valid central location determination is carried out on the data set (what is valid component, that is, theoretically belong to the same real target sample component).

  **d.** The principle of determination is to eliminate the center position caused by false noise at both ends of the order. Only when the center position difference between two adjacent components at the first end is less than 1/2 times F and the center position difference between two adjacent components at the tail end is less than 1/2 times F, it is considered that the remaining center positions in this order belong to the same real component.

  **e.** The average value of the reasonable set of central positions and their corresponding mean value of half-height and width F are used as the initial reference value of missing central positions and missing F caused by missing and false detection phenomena to participate in further LM optimization algorithm. The above is the screening criterion of RSE's Rclonte algorithm.

In the process of experimental testing, we have been thinking about a question, that is, is there a more simple and convenient way to deal with it, because the more simple, the more practical. Complex means difficult to understand and not easy to use. We wonder that is it possible to replace most of the a-e steps above by directly taking the median value of the same order as the reference value after descending the center position? After several experiments, we found that it was possible to take the median value directly after sorting, and it was completely possible on the data set collected in our experiment (see the end of the first paragraph on the right side from the seventh row from the left column on page 5 of the article in TGRS). However, in order to avoid data reuse, we simply measured two batches of data again. That is, the test of data set one stone in front, blade in back and data set two in TGRS.


### 2. **Dataset**: 
Experimental data part, the two are completely different, there is no overlap, the following is the detailed measurement situation. In the RSE experiment, we measured a batch of ASD reflectance data, including green leaves and yellow leaves, when doing biochemical component inversion. In addition, we added reflectance of dry soil and stones in the disk when doing waveform processing, which was used to establish a simulation data set. The measured data set is similar to and studied by that of the ISPRS article of Song Shalei's research group. We have done a large number of (blading - stone) distance segmentation experiments to test the effectiveness and robustness of the algorithm. In TGRS, we proposed the Rclonte-M algorithm, which was also tested on the above measured data set with good results. But at that time, we wondered whether the algorithm would be affected if the two targets changed positions. If not, it just shows that the algorithm is still feasible, so we measured the test data of the stone in front and the blade in back as data set I. The data set is not limited to this data scene, we add three target scenes to test our algorithm (stone - leaf - whiteboard in turn) to test our algorithm. The simulation dataset also uses a batch of data published by the USGS Center.


### 3. **Analysis**: 
In terms of experimental results, we have interpreted the results obtained by the above two algorithms on different test data sets in detail, such as the distance error graph not considered in RSE, which should be very interesting and can be intuitively displayed in TGRS.

Finally, in fact, in the discussion at the end of the penultimate paragraph on the right side of page 15 of RSE, we also specifically mentioned the differences between the two algorithms, mentioned the original intention of TGRS, and made a quotation explanation.


More citations to these two papers are vary appreciated. Thanks again!



