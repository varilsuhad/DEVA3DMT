% © 2020–2025 Deniz Varılsüha — Non-commercial research use only. See LICENSE.
% Contact: deniz.varilsuha@itu.edu.tr
%%% DEVA3DMT MAIN FILE TO RUN %%%

clear all;clc;close all;

inputmatrix='testinput.mat';  %%% INPUT MATRIX FROM COMING FROM THE GUI

outputfolder='C:\Users\denizv\Desktop\SVALBARD'; %%% OUTPUT FOLDER
outputname='firstinv'; %%% FILENAME TO WRITE EACH INVERSION STEP

DEVA3DMT_main(inputmatrix,outputfolder,outputname); %%% THE MAIN FUNCTION

