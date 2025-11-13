function [c,num_clust, mat] = signal_FINCH(data,initial_rank)
min_sim=inf;
[Affinity_,  orig_dist, ~]= clustRank(data,initial_rank);
initial_rank=[];
[Group_] = get_clust(Affinity_, [],inf);
[c,num_clust, mat]=get_merge([],Group_,data);
end