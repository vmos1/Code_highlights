# Introduction
The main goal is to use Generative models to produce 3D images of 
Here, we develop conditional GANs to produce images of size 128^3, conditioned on the sigma parameter.
# Plots

We develop a simple GAN trained on 3D images of size 64^3.

### Image comparison
Below are 2D snapshots of a set of 3D images. The input images are to the left and the GAN generated images are to the right.
3D GAN: Input images | 3D GAN: Generated images |
:-------------:|:---------------:
![2D slices of input images](https://github.com/vmos1/Code_highlights/blob/main/3_cond_GANs_cosmology/images/cgan_reference_2dslices.png)| ![2D slices of generated images](https://github.com/vmos1/Code_highlights/blob/main/3_cond_GANs_cosmology/images/cgan_generated_2dslices.png)

### Metric comparison
To get a better idea of the quality of generated images, we compare the pixel intensity and power spectrum of the images with the reference images.
3D cGAN: Pixel intensity | 3D cGAN: Power spectrum  |
:-------------:|:---------------:
![Pixel intensity](https://github.com/vmos1/Code_highlights/blob/main/3_cond_GANs_cosmology/images/cgan_pixel_hist.png) |![Power spectrum](https://github.com/vmos1/Code_highlights/blob/main/3_cond_GANs_cosmology/images/cgan_spec_rel.png)

As see above, the match between generated and reference images drops for larger pixel values and large k.

### Repository information
The Table below describes the important codes and their locations

| Name | Description |
| --- | ---|
| [cGAN_3d_pytorch/code/main.py](https://github.com/vmos1/Code_highlights/blob/main/3_cond_GANs_cosmology/cGAN_3d_pytorch/code/main.py) | main training code |
|[cGAN_3d_pytorch/analysis/1_cgan3d_analyze_pytorch_without_precompute.ipynb](https://github.com/vmos1/Code_highlights/blob/main/3_cond_GANs_cosmology/cGAN_3d_pytorch/analysis/1_cgan3d_analyze_pytorch_without_precompute.ipynb) | Notebook to analyze GAN results and view best epoch-steps |

## Summary: 
Although the results show a fair match, the above CGAN is quite unstable and hence difficult to train. We believe the main reason for this is the memory constraint of the 3D images. This restricts the largest batch size possible to 8. To improve the image quality and stability of the GAN, we need to train with larger batch sizes, which requires the implementation of model parallelism. We are currently working in this direction, using the LBANN framework.
