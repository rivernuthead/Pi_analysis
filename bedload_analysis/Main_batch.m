%% Simple script that Run the sequence of run name in Runs array
%ver 2.4 (11-Dec-2016)
%Author: Marco Redolfi
%Created on: 11-Dec
%Modified on: 01-Jul-2016

clear all
close all
clc

addpath('Subroutines')

%% Main parameters
% Runs = {'q07rgm', 'q07r1', 'q07r2', 'q07r3', 'q07r4', 'q07r5', 'q07r6', 'q07r7', 'q07r8', 'q07r9', 'q15r1', 'q15r2', 'q15r3', 'q15r4', 'q15r5', 'q15r6', 'q15r7', 'q15r8', 'q15r9', 'q15rgm2', 'q10rgm1', 'q10rgm2', 'q10r1', 'q10r2', 'q10r3', 'q10r4', 'q10r5', 'q10r6', 'q10r7', 'q10r8', 'q10r9', 'q20rgm2', 'q20r1', 'q20r2', 'q20r3', 'q20r4', 'q20r5', 'q20r6', 'q20r7', 'q20r8', 'q20r9'};
% Runs = {'q07r1', 'q07r2', 'q07r3', 'q07r4', 'q07r5', 'q07r6', 'q07r7', 'q07r8', 'q07r9', 'q15r1', 'q15r2', 'q15r3', 'q15r4', 'q15r5', 'q15r6', 'q15r7', 'q15r8', 'q15r9', 'q10r1', 'q10r2', 'q10r3', 'q10r4', 'q10r5', 'q10r6', 'q10r7', 'q10r8', 'q10r9', 'q20r1', 'q20r2', 'q20r3', 'q20r4', 'q20r5', 'q20r6', 'q20r7', 'q20r8', 'q20r9'};
Runs = {'q20r1', 'q20r2', 'q20r3', 'q20r4', 'q20r5', 'q20r6', 'q20r7', 'q20r8', 'q20r9'};
% Runs = {'q07r4', 'q07r5'};
%Runs = {'q07r4', 'q07r5'};
for i = 1:length(Runs)
  var_run = Runs{i}
  Main_batch_function(var_run)
end
