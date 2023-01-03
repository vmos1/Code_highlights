# Introduction
The main goal is to use Generative models to produce 3D images of 
Here, we develop conditional GANs to produce images of size 128^3, conditioned on the sigma parameter.
# Plots

We develop a simple GAN trained on 3D images of size 64^3.
### Metric comparison
3D GAN: Pixel intensity | 3D GAN: Power spectrum |
:-------------:|:---------------:
![Pixel intensity](https://github.com/vmos1/cosmogan_pytorch/blob/master/images/3d_hist_best.png)| ![Power spectrum](https://github.com/vmos1/cosmogan_pytorch/blob/master/images/3d_spec_best.png)


## 2D conditional GAN results
We develop a conditional GAN trained on 2D images of size $128^2$ for 3 different values of the cosmological paramter sigma.
### Metric comparison
2D cGAN: Pixel intensity | 2D cGAN: Power spectrum  |
:-------------:|:---------------:
![Pixel intensity](https://github.com/vmos1/cosmogan_pytorch/blob/master/images/2d_cgan_hist_best.png) |![Power spectrum](https://github.com/vmos1/cosmogan_pytorch/blob/master/images/2d_cgan_spec_best.png)

# Repository information
There separate directories to run GAN and cGAN in for both 2D and 3D datasets.
The Table below describes the content of various folders in this repository.

| Description | Location |
| --- | ---|
| 3d CGAN - code | code/5_3d_cgan/1_main_code |

There are jupyter notebooks to build launch scripts to run the code on cori GPUs at NERSC and GPUs on SUMMIT and to perform post-run computation of metrics for different stored images. Each folder contains a jupyter notebook to quickly test the code, a folder with the full code, and a folder with analysis codes to inspect the performance of the code. Below is an example for the 2D GAN:
| Name | Description |
| --- | ---|
| [1_basic_GAN/cosmogan_train.ipynb](https://github.com/vmos1/cosmogan_pytorch/blob/master/code/4_basic_3d_GAN/1_main_code/train_3dgan.ipynb) | Jupyter notebook to test conditional GAN |
| [1_basic_GAN/1_main_code](https://github.com/vmos1/Code_highlights/blob/main/3_cond_GANs_cosmology/cGAN_3d_pytorch/code/main.py) | main training code |
|[analysis/1_cgan3d_analyze_pytorch_without_precompute.ipynb](https://github.com/vmos1/Code_highlights/blob/main/3_cond_GANs_cosmology/cGAN_3d_pytorch/analysis/1_cgan3d_analyze_pytorch_without_precompute.ipynb) | Notebook to analyze GAN results and view best epoch-steps |
