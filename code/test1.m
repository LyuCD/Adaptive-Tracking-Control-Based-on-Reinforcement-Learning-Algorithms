function test1
p=[1 0 0 5 0 2 0 0 0 0 0 ]'

xn=3
P=reshape_p(p)
function P=reshape_p(p)
        P=zeros(xn);
        ij=0;
        for i=1:xn
            
            for j=1:i
                ij=ij+1;
                P(i,j)=p(ij);
                P(j,i)=P(i,j);
            end
        end
    end
end