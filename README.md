# HCP-local-processor
This project contains scripts that can be used to automate the download, processing, and clearing of HCP data. It was designed with two goals in mind, processing the large HCP dataset and doing it for free. This code should be modifable to run on most desktop workstations (needing ~6gb RAM), allowing users to still parse through the large quantity of HCP data without have the capability to store all the information at one given time or purchasing th dedicated harddrives from HCP. It also provides a free alternative to using pay-for cloud computing resources, such as Amazon EC2. This can be useful if you have a project idea with low-to-no budget.

Some inportant notes: 
1. You must have an HCP account to utilize the scripts (https://db.humanconnectome.org/app/template/Login.vm;jsessionid=FD214B9FBF54DF9D05374E91CF56EC8B). 
2. It is also dependent on a free portion of Amazon Web Service (AWS), follow these instructions provided by HCP in order to get the AWS software: https://wiki.humanconnectome.org/display/PublicData/How+To+Connect+to+Connectome+Data+via+AWS
3. This C code is dependent on a header file "nifti1.h" written by Bob Cox, not in any affiliation with this project. This is provided without any warranty or liability. Please see for details : http://paulbourke.net/dataformats/nii/nifti1.h
4. A modified version of the AAL atlas is included.

    Tzourio-Mazoyer, N. et al. Automated anatomical labeling of activations in SPM using a macroscopic anatomical parcellation of   the MNI MRI single-subject brain. Neuroimage 15, 273-289, doi:10.1006/nimg.2001.0978 (2002).
    
    Hermundstad, A. M. et al. Structural foundations of resting-state and task-based functional connectivity in the human brain. Proc  Natl Acad Sci U S A 110, 6169-6174, doi:10.1073/pnas.1219562110 (2013).
