close all
[X,Y] = meshgrid(-5: 0.05: 5); 
Z1 = X.*sin((Y -3)) + Y.*cos((X-1)/2); 
figure(1)
s1 = surf(X,Y,Z1);
colormap('winter')
zticks([]); xticks([]); yticks([])
s1.EdgeColor = 'none' ;
s1.EdgeAlpha = 0.5;

figure(2)
Z2 = peaks(X,Y);
s2 = surf(X,Y,Z2)
colormap('bone')
zticks([]); xticks([]); yticks([])
s2.EdgeColor = 'none' ;
s2.EdgeAlpha = 0.5;
