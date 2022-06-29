function [newX, sx, mx] = centernormalize(X,mode)
% Centering and Scaling the matrix X
[L,N]=size(X);

if mode==1        
        % Scaling ...
        sx=std(X,[],1);
        newX=X./sx(ones(L,1),:);

        % Centering ...
         mx=mean(newX,1);
         mmx=repmat(mx,[L 1]);
         newX=newX-mmx;
 
elseif mode==0
        newX=X;
        mx=mean(newX,1);
        mmx=repmat(mx,[L 1]);
        newX=newX-mmx;

elseif mode==2
        % Scaling ...
        sx=std(X,[],1);
        newX=X./sx(ones(L,1),:);
end
    

