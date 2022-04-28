%% Tutorial to make a simple casadi function
% This script explains how to make a simple casadi function version for the
% forward kinematics of a robot
clear;clc;close all;
%% imports
addpath(genpath(('../casadi')));
addpath(genpath(('../Robots'))); % function to get robot model
addpath(genpath(('../FloatingBaseDynamics3D'))); % this is where the forward kinematics function is located
addpath(genpath('../Utilities'))
import casadi.*
%% get a robot model
model = getMiniCheetahModel();
%% Create casadi symbolic variable to input
q = casadi.MX.sym('q_sym',model.NB,1);
%% Evaluate the matlab function symbolically
pf = casadi.MX.sym('pf_sym',3*model.N_GC,1);
pf = get_gc_position( model, q, [1:4]); 
%% make the casadi function
f_fwd_kin = casadi.Function('f_fwd_kin',{q},{pf});

%% Many ways to evaluate the casadi function numerically
% here are a few exaple
q = ones(model.NB,1);
num_val1 = f_fwd_kin(q); % Good if only a few inputs
num_val_cell = f_fwd_kin.call({q}); % Preferred way the casadi function has a lot of inputs/ouputs

%% Etra : print the SX algorithm of the forward kinematics
access_SX_algorithm(f_fwd_kin)