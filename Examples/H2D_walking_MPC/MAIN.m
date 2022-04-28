% basic LMPC for biped using single rigid body
% features: use casadi for problem formulation; solve with qpSWIFT
% Author: Yanran Ding
% Last modified: 2022/03

clear all;clc; close all;clear mex;
import casadi.*

addpath(genpath('../fcns_common'));
tic

%% initialization
LMPC = LMPC_H2D_Class;

p = LMPC.init_MPC();

p.playBackSpeed = 1;
p.SimDuration = 2;

%% MAIN Loop

p = LMPC.MAIN_LOOP(p);

%% Animation
animateRobot_SRB(p);
