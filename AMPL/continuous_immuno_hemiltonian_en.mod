##################################
#
# Sets I
#
##################################

set SIGMA; 	#The alphabet
set A;		#Allele set

##################################
#
# Parameters I
#
##################################

param N, integer, > 0;					#Length of the sequence
param eN, integer, > 0 < N;				#Length of epitopes
param k, integer, >= 0;					#Number of allowed mutations			
param pssm_thresh {A};					#Binding threshold threshold
param p {A}, >= 0.0, <= 1.0;			#Allele probability
param pssm {A, SIGMA, 1..eN};			#PSSMs for MHC-II 9mer epitope prediction

##################################
#
# Sets II
#
##################################

set Eij within {1..N,1..N}; 		#Indices which have to be considdered for the coupling optimization
set E within {1..N};				#Position of flexible positions with enougth sequence information
set M{1..N} within SIGMA;			#Position encoding -> all allowed mutation per position or only wt AA
set WT {1..N} within SIGMA;			#The wt sequence in the epitope regions

##################################
#
# Parameters II
#
##################################

param eij {(i,j) in Eij, M[i], M[j]};		#Coupling covarianz from EC (Evfold)
param h {i in E, SIGMA};					#Singelton energies

##################################
#
# Sets II
#
##################################

set Eij_neg := {
	(i,j) in Eij:
	(forall {ai in M[i], aj in M[j]} eij[i,j,ai,aj]<=0) and 
	(exists {ai in M[i], aj in M[j]} eij[i,j,ai,aj]< 0)
	};

set Eij_pos := {
	(i,j) in Eij:
	exists {ai in M[i], aj in M[j]} eij[i,j,ai,aj] > 0
	};

##################################
#
# Variables
#
##################################

var x {i in {1..N}, M[i]}, binary, >=0;				#Encoding of possible AAs of each position
var w {(i,j) in Eij, M[i], M[j]} ,binary, >=0;		#Flow variables for couplings
var y {A,i in {1..(N - eN+1)}}, >= 0;				#Immunogenicity

##################################
#
# Objective(s)
#
##################################

#Full Hammiltonian:
minimize coupling_energy: sum{(i,j) in Eij, ai in M[i], aj in M[j]} w[i,j,ai,aj] * eij[i,j,ai,aj] + sum{ i in E, ai in M[i]} x[i,ai]*h[i,ai];

##################################
#
# Constraints
#
##################################

#Only k mutations are allowed to minimize the number of epitopes
subject to only_k_mutations:
	sum{i in {1..N}, aa in WT[i]} (1 - x[i,aa]) = k;

#Each position has to choose exactly one AA
subject to only_one_aa {i in 1..N}:
	sum{aa in M[i]} x[i,aa] = 1;

#Network flow constraints
subject to flow_i_j_neg{ (i,j) in Eij_neg, ai in M[i]:
		exists {aj in M[j]} eij[i,j,ai,aj] != 0.0}:
		sum{aj in M[j]} w[i,j,ai,aj] <= x[i,ai];
		
subject to flow_j_i_neg{ (i,j) in Eij_neg, aj in M[j]:
		exists {ai in M[i]} eij[i,j,ai,aj] != 0.0}:
		sum{ai in M[i]} w[i,j,ai,aj] <= x[j,aj];

subject to flow_i_j_pos{(i,j) in Eij_pos, ai in M[i]}:
		sum{aj in M[j]} w[i,j,ai,aj] = x[i,ai];
		
subject to flow_j_i_pos{(i,j) in Eij_pos, aj in M[j]}:
		sum{ai in M[i]} w[i,j,ai,aj] = x[j,aj];

#Max functon constraint
subject to max_pssm_score{ a in A, i in {1..(N - eN+1)}}:
		y[a,i] >= (sum {j in 0..(eN - 1),aa in M[i+j]} x[i+j,aa] * pssm[a,aa,j+1]) - pssm_thresh[a];

end;
