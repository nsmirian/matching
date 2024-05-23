clear all
%addpath('OpticsTools/')

quads = {'QI.209.B1','QD.210.B1','QI.211.B1','QI.213.B1','QI.215.B1'}; %'QD.181.B1' 'QI.204.B1','QI.204.B1','QI.205.B1','QI.206.B1','QI.209.B1',  %'QD.181.B1','QI.204.B1',
reference_point = 'OTRB.218.B1';
%quads = {'QD.210.B1','QI.211.B1','QI.213.B1','QI.215.B1'}; % 

%quads = {'QD.392.B2', 'QD.415.B2','QD.417.B2','QD.418.B2' }; %,'QD.391.B2', 'QD.425.B2', 'QD.427.B2', 'QD.431.B2', 'QD.434.B2' ,'QD.437.B2','QD.440.B2'
%reference_point = 'OTRA.446.B2';

%quads = {'Q.1453.L3','QB.1475.L3','QB.1499.L3','QE.1542.L3','QE.1578.L3','QE.1615.L3','QE.1629.L3'}; %'QI.204.B1','QI.205.B1','QI.206.B1','QI.209.B1',  %'QD.181.B1','QI.204.B1',
%reference_point = 'OTRBW.1635.L3';

%[mat_x,mat_y] = FODOMatrix('B1',700);
%periodic_solution = fminsearch(@(twi,mat) match_cond(twi,twi,CurlyM(blkdiag(mat_x{3},mat_y{3}))),[1,0,1,0]);

Kmax=ones(1,length(quads))*5;
Kmin=-Kmax;
%twiss_in = [4,-1,4,1];
mismatchtol=0.01;
betamin=0.5;
betamax=200;
betaxymin=1;
Cxmin=0.5;
Cymin=0.5;
Cxymin=1;
iter=5;
Nbest=5;
Nmc=10;
withparfor=false;

zlist = [1000,1,0,1000,0];

% zlist = [207   , 1, 10 , 50, 1; ...
%         207   , 4, 0 , 50, 1]; 

% zlist = [430  , 1 , 40, 80 , 1 ; ...
%          430  , 4 ,  1, 20 , 1 ; ...
%          415  , 1 ,  5, 25 , 1 ; ...
%          415  , 4 ,  5, 25 , 1];

LLname = {'component_list',2};
IgnoreList_NAME1 = [];
[m,optics]=readLL(LLname,[quads reference_point],IgnoreList_NAME1);
% get reference position
if ~isnumeric(reference_point)
    if ischar(reference_point)
        refpos = m(end).z_pos;
        m=m(1:end-1);
    end
else
    refpos = reference_point;
end
twiss_goal = optics(end,[3,2,6,5]);
twiss_in   = [2,-0.5,8,2];%twiss_goal+[+1,0,0,0]; %twiss_goal*1.6;
%m(7).strength_sp=0.01

if withparfor
    cc=parcluster('local');
    cc.NumWorkers=16;
    pp=parpool(cc);
end

n=500; N=4; K=cell(1,iter);
for ii=1:iter
    disp(ii)
    k_start=rand(n,length(quads));k_start=bsxfun(@plus,bsxfun(@times,k_start,Kmax-Kmin),Kmin); 
    %k_start=rand(n,length(quads))*5e-2+[m([m.active]).strength_sp];
    KminI(:,ii)=Kmin;
    KmaxI(:,ii)=Kmax;
    mismatch=zeros(1,n); K{ii}(n,length(quads))=0; twiss=[];
    for i=1:n
        [out(i),stability] = the_alternator(m,refpos,k_start(i,:),twiss_in,twiss_goal,N,zlist,Kmax,Nmc,withparfor);
        mismatch(i) = out(i).mismatch_xy;
        K{ii}(i,:) = out(i).k_new;
        twiss = cat(3,twiss,out(i).matched);
        %display(num2str(stability));
        stabi(i,ii) = stability;
        misma(i,ii) = mismatch(i);
    end
    
    K{ii}=K{ii}(mismatch<mismatchtol+1,:);
    twiss=twiss(:,:,mismatch<mismatchtol+1);
    
    kmax=max(abs([Kmin Kmax]))*length(quads);
    K2sum=sum(K{ii}.^2,2);
%     Kok=K2sum<=kmax^2;
%     betaok=squeeze(min(min(twiss(:,[1 3],:),[],1),[],2))>=betamin&squeeze(max(max(twiss(:,[1 3],:),[],1),[],2))<=betamax;
%     betxyok=squeeze(min(prod(twiss(:,[1 3],:),2),[],1))>=betaxymin;
    
    indx2 = strcmp({m.type},'QUAD');
    Kall = [m(indx2).strength_sp]; 
    indx = [m(indx2).active]==1;
    Kall = repmat(Kall,size(K{ii},1),1);
    Kall(:,indx) = K{ii};
    
    Cx=sum(squeeze(twiss(out(i).quadindx,1,:).^2)'.*Kall.^2,2); %2:N+2:end
    Cy=sum(squeeze(twiss(out(i).quadindx,3,:).^2)'.*Kall.^2,2);
%     Cxok=Cx>=Cxmin;
%     Cyok=Cy>=Cymin;
%     
%     allok=Kok&betaok&betxyok&Cxok&Cyok;
%     K{ii}=K{ii}(allok,:);
%     twiss=twiss(:,:,allok);
    twiss_store{ii} = twiss;
    
%     K2sum=K2sum(allok)';
%     Cx=Cx(allok)';
%     Cy=Cy(allok)';
    mmis = min(mismatch)
    1+mismatchtol
    if mmis>1+mismatchtol
        disp('no match found');
        return
    else
    [~,icxy]=sort(Cx.*Cy);
    Cxy{ii}=Cx.*Cy;
    %Nindx = min([Nbest length(icxy)]);
    Kmin=min(K{ii}(icxy(1:Nbest),:));
    Kmax=max(K{ii}(icxy(1:Nbest),:));
    
    best(ii,:) = icxy(1:Nbest);
    end
    %figure(1); plot3(K{ii}(:,1),K{ii}(:,2),K{ii}(:,3),'.'); grid on; hold on; pause(0.01);
end
if withparfor, delete(pp); end
Nv=sum(mismatch<mismatchtol+1);

Kok=K2sum<=kmax^2;
betaok=squeeze(min(min(twiss(:,[1 3],:),[],1),[],2))>=betamin&squeeze(max(max(twiss(:,[1 3],:),[],1),[],2))<=betamax;
betxyok=squeeze(min(prod(twiss(:,[1 3],:),2),[],1))>=betaxymin;
Cxok=Cx>=Cxmin;
Cyok=Cy>=Cymin;

allok=Kok&betaok&betxyok&Cxok&Cyok;
K{ii}=K{ii}(allok,:);
twiss=twiss(:,:,allok);
twiss_store{ii} = twiss;

valid=1:Nv;
valid=valid(allok);
bb=intersect(best(ii,:),valid,'stable');
Nbest=length(bb);
best=best(:,1:Nbest);
%hold off
col = jet(iter);
% figure(101)
% for ii=1:iter
%     for j=1:length(quads)
%       for k=1:length(quads)  
%     subplot(length(quads),length(quads),j+(k-1)*length(quads))
%           plot(K{ii}(:,j),K{ii}(:,k),'.','color',col(ii,:));
%     hold on
%     plot(K{ii}(best(ii,1),j),K{ii}(best(ii,1),k),'x','color',col(ii,:));
%       end
%     end
%     grid on
% end
% hold off
% colormap(col)
% colorbar

NN=length(quads);
figure(201)
siz = round(sqrt( (NN^2-NN)/2 ));
sizy = siz;
if siz^2<(NN^2-NN)/2
    sizx = siz+1;
else
    sizx = siz;
end
count = 1;
for j=1:NN-1
    for k = j+1:NN
        subplot(sizy,sizx,count);
        
        count=count+1;
        for ii=1:iter
            plot(K{ii}(:,j),K{ii}(:,k),'.','color',col(ii,:));
            hold on
            plot(K{ii}(best(ii,1),j),K{ii}(best(ii,1),k),'x','color',col(ii,:));
        end
        hold off
        %title(['$$k_' num2str(j) ' / k_' num2str(k) '$$'],'interpreter','latex');
        xlabel(['$$k_' num2str(j) ' \left[\textrm{m}^{-2}\right]$$'],'interpreter','latex');
        ylabel(['$$k_' num2str(k) ' \left[\textrm{m}^{-2}\right]$$'],'interpreter','latex');
        grid on
    end
end


% figure(103)
% for ii=1:iter
%     for j=1:length(quads)
%       for k=1:length(quads)  
%     subplot(length(quads),length(quads),j+(k-1)*length(quads))
% %           plot(K{ii}(:,j),K{ii}(:,k),'.','color',col(ii,:));
%           scatter(K{ii}(:,j),K{ii}(:,k),[],Cxy{ii},'filled');
%     hold on
%     plot(K{ii}(best(ii,1),j),K{ii}(best(ii,1),k),'x','color',col(ii,:));
%       end
%     end
%     grid on
% end
% hold off
% colormap(col)
% colorbar

% figure(102)
% for ii=1:iter
%     plot(repmat(out(1).z',1,Nbest),squeeze(twiss_store{ii}(:,1,best(ii,:))),'color',col(ii,:))
%     hold on
%     plot(repmat(out(1).z',1,Nbest),squeeze(twiss_store{ii}(:,3,best(ii,:))),'color',col(ii,:))
% end
% c=colorbar;
% 
% plot(out(1).z',squeeze(twiss_store{iter}(:,1,best(end,1))),'-b','linewidth',3)
% plot(out(1).z',squeeze(twiss_store{iter}(:,3,best(end,1))),'-r','linewidth',3)
% hold off

figure(202)
% pos = get(202,'position');
% set(202,'position',[100,100,pos(3)*2,pos(3)*2]);
subplot(2,1,1)
for ii=1:iter
    plot(repmat(out(1).z',1,size(twiss_store{ii},3)),squeeze(twiss_store{ii}(:,1,:)),'color',col(ii,:),'linewidth',2)
    hold on
    %plot(repmat(out(1).z',1,size(twiss_store{ii},3)),squeeze(twiss_store{ii}(:,3,:)),'color',col(ii,:))
end
hold off
ylim([0,60])
grid on
xlabel('$$\textrm{position} [\textrm{m}]$$','interpreter','latex')
ylabel('$$\beta_x [\textrm{m}/\textrm{rad}]$$','interpreter','latex')
title('horizontal beta funtion')
subplot(2,1,2)
for ii=1:iter
    %plot(repmat(out(1).z',1,size(twiss_store{ii},3)),squeeze(twiss_store{ii}(:,1,:)),'color',col(ii,:))
    
    plot(repmat(out(1).z',1,size(twiss_store{ii},3)),squeeze(twiss_store{ii}(:,3,:)),'color',col(ii,:),'linewidth',2)
   hold on
end
hold off
ylim([0,60])
grid on
xlabel('$$\textrm{position} [\textrm{m}]$$','interpreter','latex')
ylabel('$$\beta_y [\textrm{m}/\textrm{rad}]$$','interpreter','latex')
title('vertical beta funtion')
%c=colorbar;

% 
% figure(205)
% test = [K{1};K{2}];
% %plot(test(:,2),test(:,4),'.')
% opts = statset('Display','final');
% [idx,C] = kmeans(test,3,'Distance','cityblock', 'Replicates',5,'Options',opts);
% plot(test(idx==1,2),test(idx==1,4),'r.','MarkerSize',12)
% hold on
% plot(test(idx==2,2),test(idx==2,4),'b.','MarkerSize',12)
% plot(test(idx==3,2),test(idx==3,4),'c.','MarkerSize',12)
% plot(C(:,1),C(:,2),'kx',...
%      'MarkerSize',15,'LineWidth',3)
% legend('Cluster 1','Cluster 2','Cluster 3','Centroids',...
%        'Location','NW')
% title 'Cluster Assignments and Centroids'
% hold off
% 
% 

NN=length(quads);
figure(208)
siz = round(sqrt( (NN^2-NN)/2 ));
sizy = siz;
if siz^2<(NN^2-NN)/2
    sizx = siz+1;
else
    sizx = siz;
end
count = 1;
for j=1:NN-1
    for k = j+1:NN
        subplot(sizy,sizx,count);
        
        count=count+1;
%         for ii=1:iter
%             plot(K{ii}(:,j),K{ii}(:,k),'.','color',col(ii,:));
%             %hold on
%             %plot(K{ii}(best(ii,1),j),K{ii}(best(ii,1),k),'x','color',col(ii,:));
%         end
Ktemp = K{1};
cxytemp = Cxy{1};
for i=2:length(K)
   Ktemp = [Ktemp; K{i}];
   cxytemp = [cxytemp ; Cxy{i}];
end
   
   indx = cxytemp < 1e9;
   
         scatter(Ktemp(indx,j),Ktemp(indx,k),15,cxytemp(indx),'filled');
         hold on
         plot(K{end}(best(end,1),j),K{end}(best(end,1),k),'xr','markersize',10,'linewidth',2);
        hold off
        %title(['$$k_' num2str(j) ' / k_' num2str(k) '$$'],'interpreter','latex');
        xlabel(['$$k_' num2str(j) ' \left[\textrm{m}^{-2}\right]$$'],'interpreter','latex');
        ylabel(['$$k_' num2str(k) ' \left[\textrm{m}^{-2}\right]$$'],'interpreter','latex');
        grid on
    end
end

centr = mean(Ktemp);
dista = sqrt(sum((Ktemp - centr).^2,2));
figure(301)
indx = dista < 10;
hist(dista(indx),500)
figure(302)
count = 1;
for j=1:NN-1
    for k = j+1:NN
        subplot(sizy,sizx,count);
        
        count=count+1;
        %         for ii=1:iter
        %             plot(K{ii}(:,j),K{ii}(:,k),'.','color',col(ii,:));
        %             %hold on
        %             %plot(K{ii}(best(ii,1),j),K{ii}(best(ii,1),k),'x','color',col(ii,:));
        %         end
        
        scatter(Ktemp(indx,j),Ktemp(indx,k),15,dista(indx),'filled');
        hold on
        %plot(K{end}(best(end,1),j),K{end}(best(end,1),k),'xr','markersize',10,'linewidth',2);
        plot(centr(j),centr(k),'xr','markersize',10,'linewidth',2);
        hold off
        %title(['$$k_' num2str(j) ' / k_' num2str(k) '$$'],'interpreter','latex');
        xlabel(['$$k_' num2str(j) ' \left[\textrm{m}^{-2}\right]$$'],'interpreter','latex');
        ylabel(['$$k_' num2str(k) ' \left[\textrm{m}^{-2}\right]$$'],'interpreter','latex');
        grid on
    end
end

thre = [0,0.5,2,4,5,6,7,10];
coll = lines(length(thre));
figure(303)
for i=1:length(thre)-1
    indx = find((dista < thre(i+1)).*(dista>thre(i)));
    scatter(Ktemp(indx,2),Ktemp(indx,4),15,coll(i,:),'filled');
    hold on
end
hold off



figure(305)
link = clusterdata(Ktemp,'maxclust',4,'linkage','single','distance','hamming');
    
for i=1:max(link)
    indx = find(link==i);
    plot(Ktemp(indx,2),Ktemp(indx,4),'.','markersize',10);
    hold on
end
hold off


figure(505)
j=2;k=4;
col =  lines(iter);
for ii=1:iter
           %scatter(K{ii}(:,2),K{ii}(:,4),15,Cxy{ii}(:),'filled');
           plot(K{ii}(:,j),K{ii}(:,k),'.k','color',col(ii,:));
           hold on
           plot([KmaxI(j,ii),KmaxI(j,ii),KminI(j,ii),KminI(j,ii),KmaxI(j,ii)], ...
                [KminI(k,ii),KmaxI(k,ii),KmaxI(k,ii),KminI(k,ii),KminI(k,ii)],'color',col(ii,:))
end
hold off
colormap(col)
%colorbar
xlabel(['$$k_' num2str(j) ' \left[\textrm{m}^{-2}\right]$$'],'interpreter','latex');
ylabel(['$$k_' num2str(k) ' \left[\textrm{m}^{-2}\right]$$'],'interpreter','latex');

figure(506)
for ii=2:iter
           %scatter(K{ii}(:,2),K{ii}(:,4),15,Cxy{ii}(:),'filled');
           plot(K{ii}(:,j),K{ii}(:,k),'.k','color',col(ii,:));
           hold on
           plot([KmaxI(j,ii),KmaxI(j,ii),KminI(j,ii),KminI(j,ii),KmaxI(j,ii)], ...
                [KminI(k,ii),KmaxI(k,ii),KmaxI(k,ii),KminI(k,ii),KminI(k,ii)],'color',col(ii,:))
end
hold off
colormap(col)
%colorbar
xlabel(['$$k_' num2str(j) ' \left[\textrm{m}^{-2}\right]$$'],'interpreter','latex');
ylabel(['$$k_' num2str(k) ' \left[\textrm{m}^{-2}\right]$$'],'interpreter','latex');

figure(507)
for ii=4:iter
           %scatter(K{ii}(:,2),K{ii}(:,4),15,Cxy{ii}(:),'filled');
           plot(K{ii}(:,j),K{ii}(:,k),'.k','color',col(ii,:));
           hold on
           plot([KmaxI(j,ii),KmaxI(j,ii),KminI(j,ii),KminI(j,ii),KmaxI(j,ii)], ...
                [KminI(k,ii),KmaxI(k,ii),KmaxI(k,ii),KminI(k,ii),KminI(k,ii)],'color',col(ii,:))
end
hold off
colormap(col)
colorbar
xlabel(['$$k_' num2str(j) ' \left[\textrm{m}^{-2}\right]$$'],'interpreter','latex');
ylabel(['$$k_' num2str(k) ' \left[\textrm{m}^{-2}\right]$$'],'interpreter','latex');