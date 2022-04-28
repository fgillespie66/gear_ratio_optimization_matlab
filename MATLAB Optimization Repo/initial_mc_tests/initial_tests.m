clc; clear all; close all;

%first generate the model
model = getMiniCheetah2DModel()
%second generate the showmotion (are these called motion priors?)
model_sm = buildMiniCheetah2DShowMotion(model)
%do the showmotion thingy
showmotion(model_sm)