clc;clear;
load matlab1.mat
refPos = [NAV(1).posX; NAV(1).posY; NAV(1).posZ];
plotGNSSSpheres(RAW_Pos, refPos, 1, 5)