

# Planar Homography Tracking for Underwater Robotics

![imagen](https://github.com/YosefGuevara012/VSLAM-MIR-MASTER/assets/54146941/799f93a8-b284-4512-99da-2e291d79efe9)


## Overview
This project focuses on the implementation of a direct planar homography tracking algorithm designed to track planar patches across sequences of images. It is particularly useful for stabilizing an Autonomous Underwater Vehicle (AUV) in relation to its environment, facilitating tasks such as underwater inspection. The codebase is provided, but completion of some core functions for the planar Homography tracking algorithm is required.

## Objective
The main objective is to gain hands-on experience with the direct planar homography tracking algorithm, enabling the sequential determination of a planar patch's location throughout an image sequence through a non-linear iterative minimization procedure.

## Key Features
- **Patch Selection:** Selection of a patch in the first image as a reference patch.
- **Initial Homography Estimation:** Initialization with the Homography H = I (identity).
- **Error Computation:** Computing the error between the reference patch and the warped patch.
- **Homography Update:** Calculating the update of H according to the pseudo-inverse of the Jacobian times the error.
- **Stopping Criterion:** Implementation of a stopping criterion for the non-linear iterative minimization.

## License
This project is licensed under [specify the license], which allows for use, modification, and distribution under specific terms.
This format includes sections for an overview, objectives, ke
