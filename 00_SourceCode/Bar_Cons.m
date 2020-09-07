%% Constitutive relationship of bars
% This function calculate the stress and material stiffness once the
% strain and material properties are given
%
% Input: 
%       Ex: the strain vector for bars;
%       barType: type of bars;
%       panelE: the Young's modulus of panel;
%       creaseE: the Young's modulus of crease;
% Output: 
%       Sx: the strain vector of the bar elements;
%       C: tangent stiffness of the bar elements;
%


function  [Sx,C]= Bar_Cons(barType,Ex,panelE,creaseE)
    Sx=zeros(size(Ex));
    C=zeros(size(Ex));

    A=size(Ex);
    N=A(1);
    for i=1:N
        if or( barType(i)==1, barType(i)==5)
            E=panelE;
        else
            E=creaseE;
        end      
        Sx(i)=E*Ex(i);
        C(i)=E;
    end
end