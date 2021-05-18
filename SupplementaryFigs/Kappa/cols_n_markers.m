% assign colors according to the topology of first layer
% and markers according to the topology of second layer

% if inet1 ==1
%     col_one = [95 158 160]/255;
% elseif inet1 ==2
%     col_one = [186,181,147]/255;
% elseif inet1 == 3
%     col_one = [167,196,139]/255;
% else
%     col_one = [231,130,162]/255;
% end
      
if jnet2 ==1
    mk = 'o';
elseif jnet2 ==2
    mk = 'sq';
elseif jnet2 == 3
    mk = 'd';    
else
    mk = '<'; 
end


if jnet1 ==1
    col_one = [0 0.4 1];
elseif jnet1 ==2
    col_one = [0.8 0.6 0];
elseif jnet1 == 3
    col_one = [0 0.6 0.2];
else
    col_one = [0.85 0.3 0.7];
end


