function [newX, sx, mx] = centernormalize3(X)
% Centering and Scaling the 2-D matrix X
[a,b,c]=size(X);
newX=zeros(size(X));
% tempX=reshape(X,a*b,c);
  for n=1:c
        % Scaling ...
        sx=nanstd(X(:,:,n),[],1);
        newXt=X(:,:,n)./sx(ones(a,1),:);

        % Centering ...
         mx=nanmean(newXt,1);
         mmx=repmat(mx,[a 1]);
         newX(:,:,n)=newXt-mmx;

% newX(:,:,n)=reshape(newX,a,b,c);         
  end