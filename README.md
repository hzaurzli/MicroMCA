# MicroMCA

## 1. Install packages
```
devtools::install_github('hzaurzli/PhagesLDA')
library(PhagesLDA)
```
**We recommend using R 4.0+ for operation**

## 2. MCA dimensionality reduction
***Step 1: Loading data***
```
data(data)
load(data)

data(info)
load(info)
```

Species abundance matrix or genes abundance matrix

***For species abundance matrix:***
| |  sample_1  |  sample_2  |  sample_3  |  sample_4 | sample_5 | sample_6 | ... | sample_m |
|  ----  | ----  |  ----  | ----  |  ----  | ----  | ---- | ---- | ---- |
| otu_1  | 8 | 9 | 11 | 4 | 2 | 9 | ... | 5 |
| otu_2  | 5 | 9 | 9 | 2 | 3 | 3 | ... | 6 |
| otu_3  | 12 | 6 | 17 | 16 | 16 | 17 | ... | 1 |
| otu_4  | 17 | 18 | 16 | 15 | 1 | 2 | ... | 6 |
| otu_5  | 15 | 16 | 14 | 2 | 6 | 5 | ... | 2 |
| ... | ... | ... | ... | ... | ... | ... | ... | ... |
| otu_n  | 7 | 8 | 5 | 2 | 6 | 17 | ... | 15 |

Rows of example data represent different otu and each column represents a different samples, count file

## 3. run MCA
```
mca = run_mca_dim(data)

# get otu' coordinates
mca$otusCoordinates

         MCA_1        MCA_2        MCA_3       MCA_4        MCA_5       MCA_6       MCA_7
otu_1  -1.5925660 -2.857963094 -1.511118115  0.17196161 -3.181968810  0.25185437  2.30685794
otu_2  -1.5437704  1.764664302 -2.170469358 -0.79385246 -0.011978713 -3.46360031 -0.09013889
otu_3  -1.5696064  0.198178377  2.249674402 -2.97040708  1.979103165 -0.20666472  2.73900827
otu_4  -1.5895893  0.629051807 -3.371606738  1.08210121  2.574874060  2.48853279  0.55206905
otu_5  -1.5562366 -1.740823626 -0.623157505 -2.40506036  0.229408382  1.07037971 -2.09260204
otu_6  -1.6071425  2.119264517  1.692760563  2.39293084 -1.003636548  0.64962992  1.44192970
otu_7  -1.5948467  0.004500578  1.709163408  0.35658022  0.211276356  2.09727916 -1.52134247
otu_8  -1.5874755  1.463384817  0.117560155 -1.08094059 -2.276088943  0.14253968 -2.11747230
otu_9  -1.5839085 -2.951441255  1.348904311  2.54116014  1.698182098 -2.20144582 -0.79715741


# get sample' coordinates
mca$samplesCoordinates

               MCA_1        MCA_2        MCA_3        MCA_4        MCA_5        MCA_6        MCA_7
sample_1    0.9730139  0.830545531  0.320451335 -1.344590674 -0.014614621 -0.582990146  0.491429720
sample_2    0.9780434  0.327198356 -0.050982925 -0.180488766  2.467342301  0.673774248 -0.354300658
sample_3    1.0043895  0.156891301 -0.877287103  0.426352378  0.510133080  0.478078929  2.217724021
sample_4    1.0257956 -0.100394991  1.162223303  0.771166429  1.064812281  1.793341498  0.999954606
sample_5    0.9972374  0.514116568 -0.111990998 -0.450677642 -0.898035427  1.930874624 -1.358022678
sample_6    0.9993955 -0.102456802  0.511406163  0.872346395  0.438609451  1.812908809 -0.847806423
sample_7    1.0087608 -0.333682961  0.758616098  1.723396713 -0.270310759 -0.228174062  0.019884597
sample_8    0.9948649 -0.339096546  0.062332056 -0.503456126  1.462279103  0.941203582  0.620752364
sample_9    0.9947220  0.149905823 -0.104547511 -1.454359909 -1.347907657  1.541504554  0.270319849
```
Visualize the results of dimensionality reduction：

![](https://github.com/hzaurzli/MicroMCA/assets/47686371/b6a8e6f5-043f-4984-8e6e-728ce61246ee)


## 4. Obtain the representative features of each phage

```
# Calculate the representative features
## Take the first 10-dimensional feature vectors after MCA dimensionality reduction for calculation

imp = MicroMCA::calculate_importance(dat = mca,info = info)

imp
$otu1
       C        A        B 
3.277906 3.582251 3.668102 

$otu2
       A        B        C 
3.621278 3.637697 3.764766 

$otu3
       B        C        A 
2.828143 2.895644 3.221709 
```
The results in imp reflect the dominant otu in each group, with smaller values reflecting that otu is more dominant in that group
