function [R,T,AA,BB,CC,DD,EE,FF] = WA_reflection_beam_area(beta,k1,k2)


AA = [1 1; -1i*k1 -k1];
BB = [1 1; 1i*k1 k1];
CC = [1 1; -1i*k2 -k2];
DD = [1i*k1 -k1; -1 1];
EE = [-1i*k1 k1; -1 1];
FF = beta^2*[1i*k2 -k2; -1 1];

R = (EE-FF*CC^-1*BB)^-1*(FF*CC^-1*AA-DD);
T = (FF-EE*BB^-1*CC)^-1*(-EE*BB^-1*CC+DD);