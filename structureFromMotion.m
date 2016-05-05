function [R,S] = structureFromMotion()
close all
%STRUCTUREFROMMOTION Summary of this function goes here
%   Detailed explanation goes here
MMatrix = load(fullfile(pwd,'assignment4_part2_data','measurement_matrix.txt'));
epsilon = 0.00001;
t  = mean(MMatrix, 2);
WH = bsxfun(@minus,MMatrix,t);
F = size(MMatrix,1)/2;
[O1,S,O2] = svd(WH);
O2 = O2';
O1p = O1(:,1:3);
Sp = S(1:3,1:3);
O2p = O2(1:3,:);
Rh = O1p * sqrt(Sp);
Sh = sqrt(Sp) * O2p;
G = createG(Rh);
c = ones(3*F,1);
c(2*F +1 : end) = 0;
Gs = pinv(G);
I  = Gs * c;
L = [I(1:3),I([2;4;5]),I([3,5,6])];
L = (L + L') /2;
[Ul,Sl,Vl] = eig(L);
Sl(Sl<0) = epsilon;
Q = Ul * sqrt(Sl);

R = Rh * Q;
S = inv(Q) * Sh;

%S = R(1,:) * S;
%R = R * R(1,:)';
axis equal;
plot3(S(1,:),S(2,:),S(3,:),'Marker','.','MarkerFaceColor','r','LineStyle','none');
Mhat = R*S;
    dif = WH-Mhat;
    mse = sum(dif.^2,2)/size(MMatrix,2);
    msec = mse(1:F,:)+mse(F+1:end,:);
    figure;
    plot(1:size(msec,1),msec);
    title('MSE per Frame')


end
function Gi = createGi(i,j)
    if nargin == 1
        j = i;
    end
    Gi = zeros(size(i,1),6);
    Gi = [repmat(i(:,1),1,3), repmat(i(:,2),1,3)];
    jcomp = [j(:,1),j(:,2),j(:,3),j(:,2),j(:,3),j(:,3)];
    Gi = Gi .* jcomp;
    addComp = [j(:,1)*0 , i(:,2).* j(:,1),...
              i(:,3).* j(:,1), i(:,2)*0, i(:,3).* j(:,2), i(:,2).* 0];
    Gi = Gi + addComp;
    
end

function G = createG(Rh)
    f = size(Rh,1)/2;
    G = [createGi(Rh(1:f,:)); createGi(Rh(f+1:end,:)); createGi(Rh(1:f,:), Rh(f+1:end,:))];
end