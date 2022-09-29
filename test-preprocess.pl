#!/usr/bin/env swipl
%%% -*- mode : prolog -*-

:- include("elbo.pl").

observed(c).
p([c,b,a],[z,y,x]).
q([c,b,a],[z,y,x]).
p([z1],[z0,a]).
