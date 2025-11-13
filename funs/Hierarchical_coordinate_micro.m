function [Y]  = Hierarchical_coordinate_micro(Y, Z, c)
Y_num=sum(Y,1);
Z_num=size(Z{1},2);
y_kT_Z=zeros(c,Z_num);
yZZy=zeros(c,1);
for k=1:c
    y_kT_Z(k,:)=Y(:,k)'*Z{1,k};
    yZZy(k)=sum(y_kT_Z(k,:).^2);
end
[~, S] = max(Y, [], 2);
S_last = 0;
obj=[];
count=1;
while any(S_last ~= S)
    S_last = S;
    for i = 1 : size(Y, 1)
        s = S(i);
        if Y_num(s) == 1
            continue
        end
        val_i=zeros(c,1);
        for k=1:c
            temp4=Z{1,k};
            if s~=k
                r1=(yZZy(k)+2*y_kT_Z(k,:)*temp4(i,:)' ...
                    +temp4(i,:)*temp4(i,:)')/(Y_num(k)+1);
                r3=yZZy(k)/Y_num(k);
                val_i(k)=r1-r3;
            else
                r1=yZZy(k)/Y_num(k);
                r3=(yZZy(k)-2*y_kT_Z(k,:)*temp4(i,:)' ...
                    +temp4(i,:)*temp4(i,:)')/(Y_num(k)-1);
                val_i(k)=r1-r3;
            end
        end
        [~, s_ba] = max(val_i);
        if s_ba ~= s
            Y(i,:)=0;
            Y(i,s_ba)=1;
            S(i)=s_ba;
            Y_num=sum(Y,1);

            y_kT_Z(s_ba,:)=Y(:,s_ba)'*Z{1,s_ba};
            y_kT_Z(s,:)=Y(:,s)'*Z{1,s};

            yZZy(s_ba)=sum(y_kT_Z(s_ba,:).^2);
            yZZy(s)=sum(y_kT_Z(s,:).^2);
        end
    end
    count = count+1;
    if count>20
        break;
    end


    
%     obj1=0;
%     for k=1:c
%         obj1=obj1+yZZy(k)/Y_num(k);
%     end
%     obj=[obj,obj1];

end
% figure(1)
% plot(1:length(obj),obj);
end