function h=pltCuboid(l1,l2,l3,g,theta,fay)
    %----------[Color map]------------------------------------------------------------
    cm(1,:)=  [     0    0.4470    0.7410];
    cm(2,:)=  [0.8500    0.3250    0.0980];
    cm(3,:)=  [0.9290    0.6940    0.1250];

%     l1=100; l2=50; l3=10;
    vert = [
        l1 -l2 -l3;
        l1 -l2 l3;
        -l1 -l2 l3;
        -l1 -l2 -l3;

        l1  l2 -l3;
        l1  l2 l3;
        -l1 l2 l3;
        -l1 l2 -l3;        ];
    fac = [1 2 3 4;
        5 6 2 1;
        8 7 6 5;
        4 3 7 8;
        2 6 7 3;
        1 5 8 4];
    
    X=vert(:,1);Y=vert(:,2);Z=vert(:,3);
    [rX,rY,rZ]=rotatez(X,Y,Z,theta);
    X=rX;Y=rY;Z=rZ;
    [rX,rY,rZ]=rotatey(X,Y,Z,fay);
    vert(:,1)=rX;vert(:,2)=rY;vert(:,3)=rZ;
    
    vert(:,1)=vert(:,1)+g(1);
    vert(:,2)=vert(:,2)+g(2);
    vert(:,3)=vert(:,3)+g(3);
    
    h=patch('Vertices',vert,'Faces',fac,...
      'FaceVertexCData',cm(3,:),'FaceColor','flat','FaceAlpha',.1,'LineWidth',0.5);
    xlabel('x');ylabel('y');zlabel('z');    
    axis equal;
    view(3);axis vis3d;
    box on; grid on;
    
end


%-------------------------------------------------------------------------- 
function [rX,rY,rZ]=rotatey(X,Y,Z,theta)
    ry(1,1)=cos(pi.*(theta)/180.);   ry(1,2)=0.;     ry(1,3)=sin(pi.*(theta)/180.);
    ry(2,1)=0.;                      ry(2,2)=1.;     ry(2,3)=0.;
    ry(3,1)=-sin(pi.*(theta)/180.);  ry(3,2)=0.;     ry(3,3)=cos(pi.*(theta)/180.);
    rX=ry(1,1).*X+ry(1,2).*Y+ry(1,3).*Z;
    rY=ry(2,1).*X+ry(2,2).*Y+ry(2,3).*Z;
    rZ=ry(3,1).*X+ry(3,2).*Y+ry(3,3).*Z;
end

%-------------------------------------------------------------------------- 
function [rX,rY,rZ]=rotatez(X,Y,Z,fay)
    rz(1,1)=cos(pi.*(fay)/180.); rz(1,2)=-sin(pi.*(fay)/180); rz(1,3)=0.;
    rz(2,1)=sin(pi.*(fay)/180.); rz(2,2)=cos(pi.*(fay)/180);  rz(2,3)=0.;
    rz(3,1)=0.;                  rz(3,2)=0.;                  rz(3,3)=1.;
    rX=rz(1,1)*X+rz(1,2)*Y+rz(1,3)*Z;
    rY=rz(2,1)*X+rz(2,2)*Y+rz(2,3)*Z;
    rZ=rz(3,1)*X+rz(3,2)*Y+rz(3,3)*Z;
end


