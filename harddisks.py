#!/usr/bin/python3.8
# -*- coding: utf-8 -*-

# This code computes the Betti numbers of the space config(n, w)
# of configurations of n unit disks in an infinite strip of width w,
# using the critical cells of a discrete Morse function constructed in
# "Configuration spaces of disks in a strip, twisted algebras,
# persistence, and other stories" by Fedor Manin and Hannah Alpert

# Enjoy! --Fedor Manin

import itertools
import math

# weights must be sorted for this to work
def bettis_of_weighted_no_k_equal_space(n, w, weights):
    bettis=[0]*n
    if weights[n-1]>w:
        return bettis
    # given a permutation, its grouping into blocks of a critical cell is unambiguous! (if we know w)
    for perm in itertools.permutations(range(n)):
        d=is_critical(perm, n, w, weights)
        if d != -1:
            #print("Critical cell "+str(perm)+" has dimension "+str(d))
            bettis[d]+=1
    return bettis

# this commented code computes all permutations as a giant list,
# obviously using itertools is better

#def all_permutations(n):
#    return all_permutations_starting_with([], list(range(n)))

#def all_permutations_starting_with(prefix, remaining):
#    if len(remaining)==0:
#        return [prefix]
#    perms=[]
#    for i in range(len(remaining)):
#        perms=perms+all_permutations_starting_with(prefix+[remaining[i]],remaining[:i]+remaining[i+1:])
#    return perms

# return -1 if not critical, dimension of cell otherwise
def is_critical(perm, n, w, weights):
    last=n
    cur_weight=0
    size_of_cur_filter=0
    d=0
    for cur in perm:
        # do we have to wrap up the current filter because it's reached max weight?
        if cur_weight>w:
            if cur_weight-min_weight>w:
                return -1
            d+=size_of_cur_filter-1
            cur_weight=0
            min_weight=weights[cur]
            size_of_cur_filter=0
        elif cur<last:
            # ...or are we ending a filter prematurely?
            if size_of_cur_filter !=0:
                return -1
            # ...or are we just riding along with singletons?
            else:
                cur_weight=0
                min_weight=weights[cur]
        # ...or are we adding on to the current filter?
        else:
            size_of_cur_filter+=1
        cur_weight+=weights[cur]
        last=cur
    if size_of_cur_filter != 0:
        if cur_weight<=w or cur_weight-min_weight>w: return -1
        else: d+=size_of_cur_filter-1
    return d

# weights must be sorted for this to work
# output is triply nested list: bettis[dimension][start of bar][end of bar]
# in output, "n+1" should be interpreted as infinity
def persistence_bars_of_weighted_no_k_equal_space(n, weights):
    big_n=sum(weights)
    bettis=[[[0 for i in range(big_n+2)] for j in range(big_n+1)] for k in range(n)]
    for perm in itertools.permutations(range(n)):
        intervals=critical_cell_bars(perm, n, weights)
        for (interval, d) in intervals:
            #print("Critical cell "+str(perm)+" has dimension "+str(d))
            bettis[d][interval[0]][interval[1]]+=1
    return bettis

# Return a list of persistence bars corresponding to critical cells
# which map to the permutation "perm" when you forget about the bars.
# For a given width w, there is at most one such cell, so we compute it
# by iterating over w.
def critical_cell_bars(perm, n, weights):
    intervals=[]
    w=weights[-1] # the least width where you get something nontrivial
    big_n=sum(weights)
    while w<big_n+1:
        last=n # previous point -- initialize it to be really big
        cur_weight=0
        size_of_cur_filter=0
        d=0 # dimension of the cell if it's critical, -1 otherwise
        next_w=big_n+1
        for cur in perm:
            # do we have to wrap up the current filter?
            if cur_weight>w:
                if cur_weight-min_weight>w: # the filter is wider than the strip!
                    d=-1
                    break
                if next_w>cur_weight:
                    next_w=cur_weight # when this filter falls apart
                d+=size_of_cur_filter-1
                cur_weight=0
                min_weight=weights[cur]
                size_of_cur_filter=0
            elif cur<last:
                # ...or are we ending a filter prematurely?
                if size_of_cur_filter !=0:
                    d=-1
                    break
                # ...or are we in the middle of some singletons?
                else:
                    cur_weight=0
                    min_weight=weights[cur]
            # ...or are we adding on to the current filter?
            else:
                size_of_cur_filter+=1
            cur_weight+=weights[cur]
            last=cur
        if d!=-1 and size_of_cur_filter != 0:
            if cur_weight<=w or cur_weight-min_weight>w:
                d=-1
            else:
                d+=size_of_cur_filter-1
                next_w=cur_weight
        if d==-1: # the cell wasn't critical
            w=w+1
        else: # we have a bar!
            intervals.append(((w, next_w), d))
            w=next_w
    return intervals

# Returns a list of tuples (p, m) where
#  * p is a partition of n into summands <= w
#  * m is the number of distinct ways of partitioning the set {1, ..., n}
#     into subsets of those sizes
def all_partitions_of_at_most_with_multiplicity(n,w):
    if w==1: return [([1]*n,1)]
    if n==0: return [([],1)]
    partitions=[]
    for i in range(n//w+1):
        sub_partitions=all_partitions_of_at_most_with_multiplicity(n-w*i,w-1)
        factor=math.prod([math.comb(n-w*j,w) for j in range(i)])//math.factorial(i)*(math.factorial(w-1)**i)
        partitions=partitions+[(sub_part+[w]*i,factor*m) for (sub_part,m) in sub_partitions]
    return partitions

def all_partitions_of_at_most(n,w):
    if w==1: return [[1]*n]
    if n==0: return [[]]
    partitions=[]
    for i in range(n//w+1):
        sub_partitions=all_partitions_of_at_most(n-w*i,w-1)
        partitions=partitions+[sub_part+[w]*i for sub_part in sub_partitions]
    return partitions

# compute the Betti numbers of config(n, w)
def bettis_of_config(n, w):
    all_weights=all_partitions_of_at_most_with_multiplicity(n,w)
    total_bettis=[0]*n
    # compute the homology of the no-k-equal space for every partition
    for weights, mult in all_weights:
        no_w_bettis=bettis_of_weighted_no_k_equal_space(len(weights), w, weights)
        dim_diff=sum(weights)-len(weights)
        for k in range(len(no_w_bettis)):
            total_bettis[k+dim_diff]+=no_w_bettis[k]*mult
    return total_bettis

# compute the Betti numbers of config(n, w) for every w
def bettis_of_each_config(n):
    bettis=[]
    for i in range(1, n+1):
        bettis.append(bettis_of_config(n,i))
    return bettis

# compute the persistence diagram of config(n, *)
def persistence_bars_of_config(n):
    all_weights=all_partitions_of_at_most_with_multiplicity(n,n)
    total_bettis=[[[0 for k in range(n+2)] for i in range(n+1)] for i in range(n)]
    # compute the persistent homology of the no-k-equal space for every partition
    for weights, mult in all_weights:
        no_w_bettis=persistence_bars_of_weighted_no_k_equal_space(len(weights), weights)
        #print("Done with weights "+str(weights)+" with multiplicity "+str(mult))
        dim_diff=sum(weights)-len(weights)
        for k in range(len(no_w_bettis)):
            for i in range(len(no_w_bettis[0])):
                for j in range(len(no_w_bettis[0][0])):
                    total_bettis[k+dim_diff][i][j]+=no_w_bettis[k][i][j]*mult
    return total_bettis

# print the persistence diagram in a nice ASCII art type format
def pretty_print_bars(n):
    bettis=persistence_bars_of_config(n)
    for k in range(n):
        print("PH bars in degree "+str(k)+":")
        for i in range(n+1):
            for j in range(n+2):
                if j==n+1:
                    j_string="===== ∞"
                else:
                    j_string=repr(j).rjust(3)+"    "
                if bettis[k][i][j] != 0:
                    print(" "*i+repr(i).rjust(2)+" "+"="*(j-i)+j_string+" "*(n+1-j)+" x"+repr(bettis[k][i][j]).rjust(n-1)+" bars")
