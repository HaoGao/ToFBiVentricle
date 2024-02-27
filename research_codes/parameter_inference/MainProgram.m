%% main program

clear all; close all; clc;

%%
'Step 1: TestAbaqusRun_preOptimization'
TestAbaqusRun_preOptimization;
%%
'Step 2: Optimization_Ca_Cb'
Optimization_Ca_Cb;
%%
'Step 3: Optimization_refine_Ca_Cb'
Optimization_refine_Ca_Cb;
%%
'Step 4: Optimization_af_bf'
Optimization_af_bf;
%%
'Step 5: Optimization_a_afs_Ca_RV'
Optimization_a_afs_Ca_RV;
