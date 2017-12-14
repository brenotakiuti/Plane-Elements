function [M,K]=EB_Beam(rho,A,L,E,I)
% M & K for the Euler-Bernoulli Bean
% Breno Ebinuma Takiuti
% 16/09/2015
% Checked 13/07/2016
% Calculate_MK.m

M = rho*A*L/420*[156    22*L     54      -13*L;
                 22*L   4*L^2    13*L    -3*L^2;
                 54     13*L     156     -22*L;
                 -13*L  -3*L^2   -22*L   4*L^2];
             
K = E*I/L^3*[12     6*L     -12     6*L;
             6*L    4*L^2   -6*L    2*L^2;
             -12    -6*L    12      -6*L;
             6*L    2*L^2   -6*L    4*L^2];