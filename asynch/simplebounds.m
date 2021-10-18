function s=simplebounds(s,Lb,Ub)
        if s(1,1) < 1 
            %s(1,1) = 1;
            s(1,1) = randi(50);
        end
        if s(1,1) > 50 
            %s(1,1)=50;
            s(1,1)=randi(50);
        end
        if s(1,2) <1 
            %s(1,2) = 1;
            s(1,2) = randi(50);
        end
        if s(1,2) > 50
            %s(1,2)=50;
            s(1,2)=randi(50);
        end  