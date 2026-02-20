function [fcz]=make_fcz(old)
% function to create z score version of 114 ROIs correlation matrix
% and use the Fisher transform to generate z score matrix
% nb: diagonal is set to zero

fc=zeros(114,114);
%make pairwise correlations
for x=1:113,
	for y=x+1:114,
		[fc(x,y),p]=corr(old(x,:)',old(y,:)');
        end;
end;

% create other half of matrix
fc_half=fc;
for x=1:113,
	for y=x+1:114,
		fc(y,x)=fc(x,y);
        end;
end;

% Fisher transform to z
fcz=.5.*log((1+fc)./(1-fc));

