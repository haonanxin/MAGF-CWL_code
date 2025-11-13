function [Y_ind,loss,Q] = main_MAGF_CWL(X,L,B,c,Y_ind)
num=size(X{1},1);
V=length(X);
D = cell(V,L,c);
dim=c;
for v=1:V
    for l=1:L
        [uu, ~, ~] = svd(B{v,l});
        B{v,l}=uu(:,1:dim);
        B{v,l}= B{v,l}./repmat(sqrt(sum(B{v,l}.^2,2)),1,dim); 
        for k=1:c
            D{v,l,k}=ones(1,dim)/sqrt(dim);
        end
    end
end

for iter=1:20
    iter
    % update Q
    f_vlk=zeros(V, L, c);
    Q=zeros(V, L, c);
    for k=1:c
        for v=1:V
            for l=1:L
                FTB=Y_ind'*B{v,l};
                FTBDk=FTB.*D{v,l,k};
                Y_ind_c=sum(Y_ind,1);
                f_vlk(v,l,:)=sum(FTB .* FTBDk,2)./Y_ind_c';
            end
        end
        Q(:,:,k)=f_vlk(:,:,k)./sqrt(sum(f_vlk(:,:,k).^2,'all'));
    end

    % update D
    for k=1:c
        for v=1:V
            for l=1:L
                D_temp = B{v,l}' * Y_ind(:, k);
                D_numerator = Q(v,l,k)*(D_temp .* D_temp);
                D_denominator = sum(Y_ind(:, k));

                D_temp2 = D_numerator / D_denominator;
                D{v,l,k} = D_temp2'./sqrt(sum(D_temp2.^2));
            end
        end
    end

    % update Y
    Z=cell(1,c);
    for k=1:c
        Bgamma=cell(V,L);
        for v=1:V
            for l=1:L
                gamma=sqrt(Q(v,l,k)*D{v,l,k});
                Bgamma{v,l}=gamma.*B{v,l};
            end
        end
        Bgamma_reshaped = reshape(Bgamma', 1, []);
        Z{1,k} = cell2mat(Bgamma_reshaped);
        clear Bgamma
    end
    Y_ind = Hierarchical_coordinate_micro(Y_ind, Z, c);

    loss(iter) = obtain_obj(Y_ind, Z, c);
     if iter > 2 && (loss(iter) - loss(iter - 1)) / loss(iter - 1) < 1e-5
        break
    end
end

end


function obj1 = obtain_obj(Y, Z, c)
obj1=0;
for k=1:c
    yZ_k=Y(:,k)'*Z{1,k};
    obj1=obj1+sum(yZ_k.^2)/sum(Y(:,k));
end

end

