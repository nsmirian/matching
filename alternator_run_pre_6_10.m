clear all
%addpath('OpticsTools/')

quads = {'QD.181.B1','QI.204.B1','QD.210.B1','QI.211.B1','QI.213.B1','QI.215.B1'}; %'QI.204.B1','QI.205.B1','QI.206.B1','QI.209.B1',  %'QD.181.B1','QI.204.B1',','QI.217.B1'
reference_point = 'OTRB.218.B1';

%[mat_x,mat_y] = FODOMatrix('B1',700);
%periodic_solution = fminsearch(@(twi,mat) match_cond(twi,twi,CurlyM(blkdiag(mat_x{3},mat_y{3}))),[1,0,1,0]);

Kmax=ones(1,length(quads))*3;
Kmin=-Kmax;
twiss_in = [4,-1,4,1];
mismatchtol=0.001;
betamin=0.5;
betamax=250;
betaxymin=1;
Cxmin=0.5;
Cymin=0.5;
Cxymin=1;
iter=5;
Nbest=5;

zlist = [214   , 4, 3 , 2000, 1; ...
         212   , 1, 0 , 15  , 1]; 


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

n=40; N=10; K=cell(1,iter);


%  twiss=[];
%  zlist = [214, 4,1,2000,0];
% for i=1:500
%     k_start=randn(length(quads),1)*50;
%     out(i) = the_alternator(m,refpos,k_start,twiss_in,twiss_goal,N,zlist);
%     mismatch(i) = out(i).mismatch_xy;
%     Ktest(i,:) = out(i).k_new;
%     twiss = cat(3,twiss,out(i).matched);
%     figure(4)
%     plot3(Ktest(:,1),Ktest(:,2),Ktest(:,3),'.'); grid on; drawnow;
% end
% return

for ii=1:iter
    disp(ii)
    k_start=rand(n,length(quads));
    k_start=bsxfun(@plus,bsxfun(@times,k_start,Kmax-Kmin),Kmin); mismatch=zeros(1,n); K{ii}(n,length(quads))=0; twiss=[];
    for i=1:n
        out(i) = the_alternator(m,refpos,k_start(i,:),twiss_in,twiss_goal,N,zlist);
        mismatch(i) = out(i).mismatch_xy;
        K{ii}(i,:) = out(i).k_new;
        twiss = cat(3,twiss,out(i).matched);
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
    
    [~,icxy]=sort(Cx.*Cy);
    %Nindx = min([Nbest length(icxy)]);
    Kmin=min(K{ii}(icxy(1:Nbest),:));
    Kmax=max(K{ii}(icxy(1:Nbest),:));
    
    best(ii,:) = icxy(1:Nbest);
    %figure(1); plot3(K{ii}(:,1),K{ii}(:,2),K{ii}(:,3),'.'); grid on; hold on; pause(0.01);
end
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
figure(101)
for ii=1:iter
    plot(K{ii}(:,1),K{ii}(:,2),'.','color',col(ii,:));
    hold on
    plot(K{ii}(best(ii,1),1),K{ii}(best(ii,1),2),'x','color',col(ii,:));
    grid on
end
hold off
colormap(col)
colorbar

figure(102)
for ii=1:iter
    plot(repmat(out(1).z',1,Nbest),squeeze(twiss_store{ii}(:,1,best(ii,:))),'color',col(ii,:))
    hold on
    plot(repmat(out(1).z',1,Nbest),squeeze(twiss_store{ii}(:,3,best(ii,:))),'color',col(ii,:))
end
plot(out(1).z',squeeze(twiss_store{iter}(:,1,best(end,1))),'-b','linewidth',3)
plot(out(1).z',squeeze(twiss_store{iter}(:,3,best(end,1))),'-r','linewidth',3)
hold off