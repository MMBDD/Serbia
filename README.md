# Serbia

1) OutlierDetection.R

  
There are typically two ways to proceed when analyzing spectra extracted from hyperspectral images. One is to gather all the spectra (pixels) corresponding to a region of interest, in this case a signal. The other is to directly make the predictions at the level of each pixel. The first option allows to drastically reduce the dimensionality of the data at  spectra level, or rows, which consequently reduces the analysis time. In this work, this first option was chosen, that is to say, all the pixels coming from the same region of interest were considered as a single spectrum.
It is important to mention, however, that not all pixels found in that region provide information, some pixels are considered as “noisy” or “bad pixels”. According to Chu Zhang et al, in 2020, “Noises care mainly caused by the environment and the sensors” https://www.sciencedirect.com/science/article/pii/S016974391930718X

According to James Burger in 2009, bad pixels show a substantial altered spectra compared to  their neighbors when detected by specific algorithms  https://journals.sagepub.com/doi/epdf/10.1255/nirn.1110. This author also provides a classification  of those abnormal pixels: “dead pixels” , “hot pixels” and “stuck pixels”. Dead pixels do not react to light.  According to Bingkai  Liu in 2023,  “Hot pixels refer to the anomalous pixel with high dark current compared to the normal pixels with moderate dark current increase after irradiation”, https://www.mdpi.com/1424-8220/23/13/6159. Finally, stuck pixels show an almost steady  intermediate value. Moreover, it is interesting to note, that some pixels are always noisy, while others only sporadically; some may show a “non- linear response to light intensity” https://journals.sagepub.com/doi/epdf/10.1255/nirn.1110, while some others behave randomly. 

In any case, these abnormal pixels behave differently from the rest and can be detected by algorithms in an unsupervised manner. In this work, the following methods for detecting these pixels were compared in relation to the results obtained in the final predictions: (i) Principal Component Analysis (PCA), (ii) Isolation Forest (IF),  (iii) PCA + IF (iv) one-class Soft Independent Modeling by Class Analogy (SIMCA), (v) PCDIST algorithm,  (vi) one-class Support Vector Machine and (vii) an univariate approach. 

Additionally, three ways of combining relevant spectra belonging to the same area of interest were studied: (i) by taking their average, (ii) by calculating their standard deviation, adding this value as an extra variable, and finally computing the average at the row level and (iii) by pretreating these samples with the Standard Normal Variate pretreatment (SNV), and then calculating the average.  In total, 21 predictions were compared to have a preliminary idea of their performance in the problem at hand. 


The results will be compared in terms of the different methods used, and also in terms of different ways to combine relevant pixels.

With this code, you will be able to transform a Spectra Matix that corresponds to one Sepal, into only One Spectrum per Sepal. A final matrix with One Spectrum per Sepal; will be used as input for the next R code (GlobalModel.R)

2) GlobalModel.R


The data sets of varieties Brioso, Cappricia and Provine were uploaded to the R code “GlobalModel.R”.

The following operations were carried out in this code:

2.1. Exploration by PCA in each data set. The outliers removed can be seen below: 
Cultivar	Outliers removed
Brioso	"2599", "4038", "983", "3222", "3955", "3651", "440", "5160"
Cappricia	"i2T11S5", "i2T5S3", "i1T1S3", "i1T1S2"
Provine	"7790", "9639", "0", "7232","10340", "1002", "4117"

2.2. 	Samples were distributed in two classes according to visual scoring. The data sets were divided into calibration (70%) and validation (30%) sets, in a representative way for each class. In other words, Class 1 was split in a proportion of 70/30 and Class 2 was also split in the same proportion 70/30

2.3 In the present work an iterative process was used to select a sparse subset of important variables before their use by the classification models. Important variables were optimized for each pretreatment, labelling and cultivar. 
The selected variables as given by CovSel were then used as input of PLSDA. The training set was split again into Validation and Tuning and the optimal number of latent variables was selected according to the lower error shown in the PLSDA model. The discriminant model calibrated on the training subset was applied to independent samples in the test set. The results were expressed with classification metrics such as sensitivity, specificity, precision, accuracy and balanced accuracy. Besides raw data, several preprocessing steps were performed and compared. Models were built in the training set using 5 to 39 selected variables by CovSel. PLSDA latent variables were optimized as well, by cross-validation on each tomato.

The iterative process was carried out in the Training Set, and can be summarized as follows: 

3.1: First of all, a pretreatment was chosen in order to remove noise from the original spectra. The following algorithms were compared: Detrend grades 1 and 2; Savitzky–Golay first and second derivatives, second polynomial degree and 9, 11, 15, 17 smoothing windows; Standard Normal Variate (SNV); and combinations of these. 

3.2: Then, the CovSel algorithm was applied to the pretreated data, 5 important variables were chosen. 

3.3: The next step consisted of dividing the Training Set again (which now contains only the important variables); in Validation and Tuning, randomly.

3.4: The selected variables as given by CovSel were then used as input of PLSDA. The Validation set was used to train the PLSDA classification model, and the Tuning set was used to validate it. This last process was carried out with the objective of choosing the optimal number of latent variables for the PLSDA algorithm. The number of latent variables that showed the least error was chosen.

3.5: Steps 3.2 to 3.4 were repeated with an increasing number of variables from 5 to 39.

3.6: Steps 3.1 to 3.5 were repeated by choosing a different pretreatment.

3.7: Steps 3.1 to 3.6 were repeated using Labeling scenarios 2 and 3. 

3.8: Steps 3.1 to 3.7 were repeated using different varieties "Cappricia", "Provine" and "Brioso".

3.9: Varieties were added together (Cappricia + Provine, Cappricia + Provine + Brioso, Cappricia + Brioso, Brioso + Provine), then steps 3.1 to 3.7 were repeated.

3.10. Steps 3.1 to 3.7 were replicated, choosing the Cappricia variety as the Training set and the Brioso variety as the Test Set.

3.11. Item 3.10 was reiterated, but choosing Provine as a Test Set.

Out of all the results obtained, the model that presented the lowest Balanced Accuracy (BA) value was chosen. 

3) PLSDA.R   Build PLSDA models with the optimized parameters using Rchemo
4) PLSDAmdatools.  Build PLSDA and SIMCA models using mdatools
5) MSC.R   Apply pretreatments to raw spectra, and plot the results.