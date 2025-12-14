clearvars
close all

C = 1.0e5*[8,2,0;2,8,0;0,0,3]/12.0;
B = {[-4,0,4,0,0,0; 0, 0, 0, -3, 0, 3; 0, -4, -3, 4, 3, 0], ...
    -[-4,0,4,0,0,0; 0, 0, 0, -3, 0, 3; 0, -4, -3, 4, 3, 0]};

u = [
    0;
    0;
    1.129111438389789e-01;
    1.963672066764852e-02;
    1.011291114383898e-01;
   -1.080019636720666e-02;
    0;
    0;
    ];

elem = [1, 2, 3;
    3, 4, 1];

numElem = size(elem,1);
dim = 2;

sigma = zeros(3,numElem);

for e=1:numElem
    rows = [dim*elem(e,1)-1, dim*elem(e,1), ...
        dim*elem(e,2)-1, dim*elem(e,2),...
        dim*elem(e,3)-1, dim*elem(e,3)];
    sigma(:,e) = C*B{e}*u(rows);
end

vonMises = sqrt(sigma(1,:).^2 + sigma(2,:).^2 ...
    - sigma(1,:).*sigma(2,:)+3*sigma(3,:).^2);

[sigma', vonMises']

  