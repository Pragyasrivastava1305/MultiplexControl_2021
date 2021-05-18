% define colors
clen = 5; 
blue = [0,0,1];
green1 = [0, 0.8, 0];
colors_p1 = [linspace(blue(1),blue(1),clen)', linspace(blue(2),green1(2),clen)', linspace(blue(3),green1(3),clen)'];

gray1 =[0.4 0.4 0.4]; 
yellow1 = [204, 153, 0]/255;  yellow2 = [255, 204, 153]/255; 
mauve1  = [204,51,153]/255;    mauve2 = [255,153,255]/255; 
green1 = [0, 153,51]/255;     green2 = [0,255,205]/255;
blue1 = [102,0,255]/255;     blue2 = [0,204,255]/255; 

% one color scheme
cvec1 = [linspace(blue(1),blue2(1),clen)', linspace(blue1(2),blue2(2),clen)', linspace(blue1(3),blue2(3),clen)'];
cvec2 = [linspace(yellow1(1),yellow2(1),clen)', linspace(yellow1(2),yellow2(2),clen)', linspace(yellow1(3),yellow2(3),clen)'];
cvec3 = [linspace(green1(1),green2(1),clen)', linspace(green1(2),green2(2),clen)', linspace(green1(3),green2(3),clen)'];
cvec4 = [linspace(mauve1(1),mauve2(1),clen)', linspace(mauve1(2),mauve2(2),clen)', linspace(mauve1(3),mauve2(3),clen)'];
% % second color scheme
% col0 = gray1; 
 col1 = yellow2; 
 col2 = green2;
 col3 = blue2; 
 col4 = mauve2; 
% 
% 
% cvec1 = [linspace(col0(1), col1(1),clen)', linspace(col0(2), col1(2),clen)', linspace(col0(3), col1(3),clen)'];
% cvec2 = [linspace(col0(1), col2(1),clen)', linspace(col0(2), col2(2),clen)', linspace(col0(3), col2(3),clen)'];
% cvec3 = [linspace(col0(1), col3(1),clen)', linspace(col0(2), col3(2),clen)', linspace(col0(3), col3(3),clen)'];
% cvec4 = [linspace(col0(1), col4(1),clen)', linspace(col0(2), col4(2),clen)', linspace(col0(3), col4(3),clen)'];
% 
 bcol1 = cvec1(3,:);  
 ycol1 = cvec2(1,:); 
 gcol1 = cvec3(1,:);  
 mcol1 = cvec4(2,:);