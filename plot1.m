%function []=plot1(r1,z1,us1,Npar,cnteta,nray)
figure(1111)
h1=axes('position',[0.1 0.1 0.8 0.8]);
box on;
hold on;
plot(gvar.rbbbs,gvar.zbbbs,'r','LineWidth',2);   %最外闭合磁面
plot(gvar.rlim,gvar.zlim,'b','LineWidth',2);    %限制器
axis equal;
xlabel('R(m)');
ylabel('Z(m)');
plot(r1(1,1),z1(1,1),'*r', 'markersize',15)
title(['f=',num2str(w11),' N_{//}=',num2str(npar)]);
LB=axis;
hold off;
h2=axes('position',[-0.15 0.1 0.75 0.8]);
% axis(h2);

set(gca,'xticklabel',[]);
set(gca,'yticklabel',[]);
hold on;
box on;
h2=plot(gvar.rbbbs,gvar.zbbbs,'r','LineWidth',2);   %最外闭合磁面
h2=plot(gvar.rlim,gvar.zlim,'b','LineWidth',2);    %限制器
axis equal;

plot(r1(1,1),z1(1,1),'*r', 'markersize',15)
axis([2.2 2.35 -0.2 0.2])
for j=1:50:nray
   
    
%     [k1,~]=min(find(r1(:,j)==0));
%     if inpolygon(r1(end,j),z1(end,j),gvar.rlim,gvar.zlim)
%         k1=length(r1(:,j));
%     else
%         k1=length(r1(:,j))-1;
    if r1(end,j)>0
        k1=length(r1(:,j))-2;
     else
         [k1,~]=min(find(r1(:,j)==0));
        k1=k1-2;
    end
    if j==1
        
        figure(1111)
        hold on;
        h1=plot(r1(1:k1,j),z1(1:k1,j),'LineWidth',2,'color','b');

        h2=plot(r1(1:k1,j),z1(1:k1,j),'LineWidth',2,'color','b');
%         figure(2)
%         hold on;
%         plot(teta1(1:k1,j),Npar1(1:k1,j),'LineWidth',2,'color','b');
%         max(teta1(1:k1,j));
%             
%        
%         figure(3)
%         hold on;
%         plot(teta1(1:k1,j),cnteta1(1:k1,j),'LineWidth',2,'color','b');
%         
        figure(2222)
        hold on;
        plot(us1(1:k1,j),Npar1(1:k1,j),'LineWidth',2,'color','b');
        figure(3333)
        hold on;
        plot(us1(1:k1,j),cnteta1(1:k1,j),'LineWidth',2,'color','b');
%         figure(6)
%         hold on;
%         plot(us1(1:k1,j),cnprim1.e(1:k1,j),'LineWidth',2,'color','b');
        

        
        
    else
        figure(1111)
        hold on;
        h1=plot(r1(1:k1,j),z1(1:k1,j),'LineWidth',2);
        h2=plot(r1(1:k1,j),z1(1:k1,j),'LineWidth',2);
        figure(2222)
        hold on;
        plot(us1(1:k1,j),Npar1(1:k1,j),'LineWidth',2);
        figure(3333)
        hold on;
        plot(us1(1:k1,j),cnteta1(1:k1,j),'LineWidth',2);
    end

    
% if j~=1
%     if k1<length(r1(:,1))-1
%         k1
% figure(1)
% hold on;
% plot(r1(1:k1,j),z1(1:k1,j),'LineWidth',2,'color','b');
% figure(2)
% hold on;
% plot(us1(1:k1,j),Npar(1:k1,j),'LineWidth',2,'color','b');
% figure(3)
% hold on;
% plot(us1(1:k1,j),cnteta1(1:k1,j),'LineWidth',2,'color','b');
%     else 
%         figure(1)
% hold on;
% plot(r1(:,j),z1(:,j),'LineWidth',2,'color','r');
% figure(2)
% hold on;
% plot(us1(:,j),Npar(:,j),'LineWidth',2,'color','r');
% figure(3)
% hold on;
% plot(us1(:,j),cnteta1(:,j),'LineWidth',2,'color','r');
%     end
% else
%     %if k1==[];
%     k2=min(length(r1(:,j))-1,k1)
%         figure(1)
%         hold on;
%         plot(r1(1:k2,j),z1(1:k2,j),'LineWidth',3,'color','k');
%         figure(2)
%         hold on;
%         plot(us1(1:k2-1,j),Npar(1:k2-1,j),'LineWidth',3,'color','k');
%         figure(3)
%         hold on;
%         plot(us1(1:k2,j),cnteta1(1:k2,j),'LineWidth',3,'color','k');  
%    % end
% end


end
figure (2222)
    box on;
    title(['f=',num2str(w11),' N_{//}=',num2str(npar)]);
    xlabel('ray tracing distance(m)');
    ylabel('N_{//}');
figure (3333)
    box on;
    title(['f=',num2str(w11),' N_{//}=',num2str(npar)]);
    xlabel('ray tracing distance(m)');
    ylabel('N_{\theta}');
    
% saveas(1111,[savePath1,'_rt.fig']);
% saveas(2222,[savePath1,'_NparS.fig']);
% saveas(3333,[savePath1,'_cntetaS.fig']);