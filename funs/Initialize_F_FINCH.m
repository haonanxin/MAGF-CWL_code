function [FF,num_view] = Initialize_F_FINCH(X,c)

% [F_set, num_cluster]=FINCH(X{1},[],[]);
% for i=2:length(num_cluster)
%     if num_cluster(i-1)>=c && num_cluster(i)<c
%         req_idx=i-1;
%     elseif num_cluster(length(num_cluster))>=c
%         req_idx = length(num_cluster) ;
%     end
% end
% dd=F_set(:,req_idx);
% F_init=req_numclust(F_set(:,req_idx),X{1},c);
% F = full(ind2vec(F_init')');

V=length(X);
Q_all=zeros(V,1);
F_all=cell(V,1);
num_view=zeros(V,1);
for v=1:V
    [F_set, num_cluster]=FINCH(X{v},[],[]);
    if num_cluster~=1
        num_view(v)=sum(num_cluster>=c);
        for i=2:length(num_cluster)
            if num_cluster(i-1)>=c && num_cluster(i)<c
                req_idx=i-1;
            elseif num_cluster(length(num_cluster))>=c
                req_idx = length(num_cluster) ;
            end
        end
        F_init=req_numclust(F_set(:,req_idx),X{v},c);
        F = full(ind2vec(F_init')');
%         A = constructW_PKN(X{v}, 10);
        A = X{v}*X{v}';
        [Q] = modularity(A,F);
        Q_all(v)=Q;
        F_all{v}=F;
    else
        Q_all(v)=-100;
        F_all{v}=0;
    end
end
[M,idx] = max(Q_all);
FF=F_all{idx};

end


function [Q] = modularity(A,S)
% A = [1,1,0;1,1,0;0,0,1]; %邻接矩阵定义，与上面的例子是一致的
% S = [1,0;1,0;0,1]; %label定义
m = sum(sum(A))/2;
k = sum(A,2);
B = A - (repmat(k,[1,size(A,1)]) .* repmat(k',[size(A,1),1])) / (2*m);
Q = 1/(2*m) .* trace(S'*B*S);
end
