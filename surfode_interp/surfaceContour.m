% 3D contour plot on surface
% Adapted from plotFS.m which is written by Tim

function surfaceContour(X,Y,Z,F,V,T,nLevel)

% figure
hold on

C=interp3(X,Y,Z,T,V(:,1),V(:,2),V(:,3));
cMin=min(C);
cMax=max(C);

for level=1:nLevel
    cLevel=cMin+(cMax-cMin)/nLevel*(level-1);
    for i=1:length(F)
        
        xA=V(F(i,1),:);
        xB=V(F(i,2),:);
        xC=V(F(i,3),:);
        
        fA=C(F(i,1));
        fB=C(F(i,2));
        fC=C(F(i,3));
        
        tmp1=(cLevel-fA)*(cLevel-fB);
        tmp2=(cLevel-fB)*(cLevel-fC);
        tmp3=(cLevel-fA)*(cLevel-fC);
        
        tmp=min([tmp1,tmp2,tmp3]);
        
        if tmp<=0
            if (tmp1<=0)&&(tmp2<=0)
                y1=(abs(fB-cLevel)*xA+abs(fA-cLevel)*xB)/...
                   (abs(fB-cLevel)+abs(fA-cLevel));
                y2=(abs(fB-cLevel)*xC+abs(fC-cLevel)*xB)/...
                   (abs(fB-cLevel)+abs(fC-cLevel));
            elseif (tmp2<=0)&&(tmp3<=0)
                y1=(abs(fC-cLevel)*xA+abs(fA-cLevel)*xC)/...
                   (abs(fC-cLevel)+abs(fA-cLevel));
                y2=(abs(fB-cLevel)*xC+abs(fC-cLevel)*xB)/...
                   (abs(fB-cLevel)+abs(fC-cLevel));
            elseif (tmp1<=0)&&(tmp3<=0)
                y1=(abs(fB-cLevel)*xA+abs(fA-cLevel)*xB)/...
                   (abs(fB-cLevel)+abs(fA-cLevel));
                y2=(abs(fA-cLevel)*xC+abs(fC-cLevel)*xA)/...
                   (abs(fA-cLevel)+abs(fC-cLevel));
            end
            p = plot3([y1(1) y2(1)],[y1(2) y2(2)],[y1(3) y2(3)],'k-',...
            'LineWidth',1,'MarkerSize',2) ;
            p.Color(4) = 0.5 ;
        end
    end
end

axis equal
axis off