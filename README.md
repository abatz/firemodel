# firemodel

System Requirements
Software required: MATLAB, Statistics Toolbox
Version tested: R2021a

Installation Instructions
Place all files in a directory 

Demo
Execute main code runcode.m 
Expected output:                 
For example, user query of 30-model median burned area for model years 72-101 (2021-2050)
>> squeeze(nanmedian(nanmedian(nanmean(burned(72:101,:,1,:),1),2),4))

ans =

   1.0133e+04

Typical runtime with files provided: 2500 seconds
