%%% -*- mode : prolog -*-

:- use_module(library(lists)).

%%% Generic utilities.

setequal(X,Y):-
    subset(X,Y),
    subset(Y,X).

union(A,B,C,R2):-
    union(A,B,R1),
    union(R1,C,R2).
union(A,B,C,D,R2):-
    union(A,B,C,R1),
    union(R1,D,R2).
union(A,B,C,D,E,R2):-
    union(A,B,C,D,R1),
    union(R1,E,R2).

rsum([X], X).
rsum([X1,X2|Xs], Sum+X1):-       % ensure more than two elements
    rsum([X2|Xs], Sum).

sum(Xs,Sum):-
    reverse(Xs,RXs),
    rsum(RXs,Sum).



%% anysubset(Set, Subset).
%% enumerate subsets of Set.
anysubset([], []).
anysubset([E|Es], [E|Sub]):-
    anysubset(Es, Sub).
anysubset([_|Es], Sub):-
    anysubset(Es, Sub).

%%% preprocessing
%% sort variables according to the names for all distributions.
%% This replaces user-defined ones with sorted ones, which requires dynamic/1.

:- dynamic p/2, q/2.
preprocess:-
    nl,
    writeln("preprocessing: reordering the variables in p"),
    findall([p(X1,Y1),p(X2,Y2)],
            (p(X1,Y1),
             sort(X1,X2),
             not(=(X1,X2)),
             sort(Y1,Y2),
             not(=(Y1,Y2))),Ps),
    foreach(member([p(X1,Y1),p(X2,Y2)], Ps),
            (retract(p(X1,Y1)),
             assertz(p(X2,Y2)),
             write("replacing "),write(p(X1,Y1)),write(" -> "),write(p(X2,Y2)),nl)),
    writeln("preprocessing: reordering the variables in q"),
    findall([q(X1,Y1),q(X2,Y2)],
            (q(X1,Y1),
             sort(X1,X2),
             not(=(X1,X2)),
             sort(Y1,Y2),
             not(=(Y1,Y2))),Ps),
    foreach(member([q(X1,Y1),q(X2,Y2)], Ps),
            (retract(q(X1,Y1)),
             assertz(q(X2,Y2)),
             write("replacing "),write(q(X1,Y1)),write(" -> "),write(q(X2,Y2)),nl)),
    compile_predicates([p/2, q/2]),
    writeln("preprocessing: done."),
    true.


%%% random variables.

:- table observables/1, variables/1, latents/1.

observables(Xs):-
    findall(X, observed(X), Xs1),
    sort(Xs1,Xs).

variables(Xs):-
    findall(X, (p(V1,V2),union(V1,V2,V),member(X,V)), Xs1),
    findall(X, (q(V1,V2),union(V1,V2,V),member(X,V)), Xs2),
    union(Xs1,Xs2,Xs3),
    sort(Xs3,Xs).

latents(Xs):-
    variables(Vars),
    observables(Obs),
    subtract(Vars,Obs,Xs).


%%% distributions

:- table generative/1, variational/1, empirical/1, split_P/4.

generative(P):-
    findall(p(Vars,Deps),p(Vars,Deps),P).
variational(Q):-
    latents(Zs),
    findall(q(Vars,Deps),(q(Vars,Deps),subset(Vars,Zs)),Q).
empirical(Q):-
    observables(Xs),
    findall(q(Vars,Deps),(q(Vars,Deps),subset(Vars,Xs)),Q).


%%% ELBO

split_P(P,Q,P1,P2,P3):-
    observables(Xs),
    latents(Zs),
    findall(p(Vars,Deps),
            (member(p(Vars,Deps),P),
             subset(Vars,Zs)),
            P0),
    findall(p(Vars,Deps),
            (member(p(Vars,Deps),P0),
             member(q(Vars,_),Q)),
            P1),
    findall(p(Vars,Deps),
            (member(p(Vars,Deps),P0),
             not(member(q(Vars,_),Q))),
            P2),
    findall(p(Vars,Deps),
            (member(p(Vars,Deps),P),
             subset(Vars,Xs)),
            P3).

elbo('E'(Sampling_Distributions, Formula), Sampling_Distributions):-
    generative(P),
    variational(AllQ),
    empirical(Emp),
    anysubset(AllQ,Q),
    split_P(P,Q,P1,P2,P3),
    union(Emp,Q,P2,Sampling_Distributions),
    %% generate the form
    findall(log(p(Vars,Deps1))/log(q(Vars,Deps2)),
            (member(p(Vars,Deps1), P1),
             member(q(Vars,Deps2), Q)),
            LogRatio),
    findall(log(p(Vars,Deps)),
            member(p(Vars,Deps), P3),
            Reconstruction),
    union(Reconstruction,LogRatio,Terms),
    sum(Terms, Formula),
    true.



%%% ensure that all distributions can be sampled

dependency_hypergraph(Distributions,HyperGraph):-
    findall(Vars-Deps, member(p(Vars,Deps),Distributions), HyperGraph_P),
    findall(Vars-Deps, member(q(Vars,Deps),Distributions), HyperGraph_Q),
    union(HyperGraph_P,HyperGraph_Q,HyperGraph).

targets(HyperGraph, Targets):-
    setof(Var, (member(Vars-_,HyperGraph),member(Var,Vars)), Targets).

:- table samplable/1, samplable/2.

samplable(Distributions):-
    dependency_hypergraph(Distributions, HyperGraph),
    targets(HyperGraph, Targets),
    samplable(HyperGraph, Samplable),
    subset(Targets, Samplable).

samplable([],[]).

samplable([Vars-[]],Vars).

samplable(HyperGraph,NewSamples):-
    member(Vars-Deps, HyperGraph),
    subtract(HyperGraph, [Vars-Deps], SubHyperGraph),
    samplable(SubHyperGraph, Samples),
    subset(Deps,Samples),
    union(Vars,Samples,NewSamples).



%% main

report_list(List):-
    forall(
        member(Elem, List),
        writeln(Elem)).


main_all_elbos:-
    nl,
    writeln("ELBOs that may including unrealizable sampling process:"),
    findall(F,elbo(F,_),Fs),
    length(Fs,N),
    write(N),writeln(" ELBOs."),
    nl,
    report_list(Fs),
    true.


main_samplable_elbos:-
    nl,
    writeln("ELBOs with a feasible sampling process:"),
    findall(F,
            (elbo(F,Distributions),
             samplable(Distributions)),
            Fs),
    length(Fs,N),
    write(N),writeln(" ELBOs."),
    nl,
    report_list(Fs),
    true.





main:-
    writeln("##################################"),
    preprocess,
    nl,
    writeln("Analyzing the generative model of:"),
    observables(Xs),
    writeln('E'(q(Xs), p(Xs))),
    nl,

    main_all_elbos,
    nl,
    main_samplable_elbos,
    nl,
    halt,
    true.

:- initialization(main).
