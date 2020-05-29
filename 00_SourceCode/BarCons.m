%% Constitutive relationship of bars
% This function calculate the stress and material stiffness once the
% strain and material properties are given

% The content in this code is based on the open-access deformable orgami 
% simulator developed by K. Liu and G. H. Paulino 
% [1] K. Liu, G. H. Paulino (2017). 'Nonlinear mechanics of non-rigid   
%     origami - An efficient computational approach.' PRSA.  

function  [Sx,C]= BarCons(BarType,Ex,PanelE,CreaseE)
    Sx=zeros(size(Ex));
    C=zeros(size(Ex));

    A=size(Ex);
    N=A(1);
    for i=1:N
        if or( BarType(i)==1, BarType(i)==5)
            E=PanelE;
        else
            E=CreaseE;
        end      
        Sx(i)=E*Ex(i);
        C(i)=E;
    end
end