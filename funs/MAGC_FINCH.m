function [B] = MAGC_FINCH(X,L,alpha)
V=length(X);
num=size(X{1},1);
P1=cell(V,L);
P2=cell(V,L);

Q=cell(V,L);
Center=cell(V,L);
for v=1:V
    Center{v,1}=X{v};
end

for l=1:L
    for v=1:V
        [idx,~, center] = signal_FINCH(Center{v,l},[]);
        Q{v,l}=full(ind2vec(idx')');
        Center{v,l+1}=center;
        if l>1
            P1{v,l}=P1{v,l-1}*Q{v,l};
            P2{v,l-1}=P1{v,l}*Q{v,l}';
        elseif l==1
            P1{v,1}=Q{v,l};
        end
    end
end

k=10;
for v=1:V
    for l=1:L-1
        m=size(Center{v,l+1},1);
        if k<m
            K1=get_knn(X{v},Center{v,l+1},k);
            K2 = P2{v,l}./sum(P2{v,l},2);
            B{v,l}=(1-alpha)*K1+alpha*K2;
        else
            K1=get_knn(X{v},Center{v,l+1},m-1);
            K2 = P2{v,l}./sum(P2{v,l},2);
            B{v,l}=(1-alpha)*K1+alpha*K2;
        end
    end
end
end






% if l==1
%     K1 = double(exp(-pdist2(X{v},Center{v,l+1}, 'euclidean').^2 / (2 * sigma^2)));
%     K2 = K1.*Q{v,1};
%     K2=K2./sum(K2,1);
%     T{v,l}=K2;
% elseif l==2
%     K1 = double(exp(-pdist2(Center{v,l}, Center{v,l+1}, 'euclidean').^2 / (2 * sigma^2)));
%     K2 = K1.*Q{v,l};
%     K2=K2./sum(K2,1);
%     T{v,l}=K2;
% else
%     K1 = double(exp(-pdist2(Center{v,l}, Center{v,l+1}, 'euclidean').^2 / (2 * sigma^2)));
%     K2 = K1.*Q{v,l};
%     K2=K2./sum(K2,1);
%     K3=T{v,l-1}*K2;
%     T{v,l}=K3./sum(K3,1);
% end
% 
% if l==1
%     R{v,l}=T{v,1};
% else
%     temp_T=T{v,l}*T{v,l}';
%     R{v,l}=T{v,1}*temp_T;
% end