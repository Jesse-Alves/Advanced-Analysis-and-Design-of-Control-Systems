clear all
close all
pkg load control

I = eye(2);

Q = I;
R = 0.1;

A1 = [0 1; -1 0];
B1 = [0; 1];

A2 = [0 1; -1 -4.5];
B2 = [0; 1];

[Klqr1, S, E] = lqr(A1, B1, Q, R);

Klqr1

[Klqr2, S, E] = lqr(A2, B2, Q, R);

Klqr2