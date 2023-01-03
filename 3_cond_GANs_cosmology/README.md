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
The Table below describes the important codes and their locations

| Name | Description |
| --- | ---|
| [cGAN_3d_pytorch/code/main.py](https://github.com/vmos1/Code_highlights/blob/main/3_cond_GANs_cosmology/cGAN_3d_pytorch/code/main.py) | main training code |
|[cGAN_3d_pytorch/analysis/1_cgan3d_analyze_pytorch_without_precompute.ipynb](https://github.com/vmos1/Code_highlights/blob/main/3_cond_GANs_cosmology/cGAN_3d_pytorch/analysis/1_cgan3d_analyze_pytorch_without_precompute.ipynb) | Notebook to analyze GAN results and view best epoch-steps |
