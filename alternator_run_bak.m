clear all
%addpath('OpticsTools/')

quads = {'QD.210.B1','QI.211.B1','QI.213.B1','QI.215.B1','QI.217.B1'}; %'QI.204.B1','QI.205.B1','QI.206.B1','QI.209.B1',
reference_point = 'OTRB.218.B1';

%[mat_x,mat_y] = FODOMatrix('B1',700);
%periodic_solution = fminsearch(@(twi,mat) match_cond(twi,twi,CurlyM(blkdiag(mat_x{3},mat_y{3}))),[1,0,1,0]);

Kmax=[2 2 2 2 2];
Kmin=-Kmax;
twiss_in = [4,-1,4,1];
mismatchtol=0.001;
betamin=0.5;
betamax=250;
betaxymin=1;
Cxmin=0.5;
Cymin=0.5;
Cxymin=1;
iter=10;
Nbest=10;

zlist = [214   , 4, 3 , 2000, 1; ...
         212   , 1, 0,15   , 1]; 


LLname = {'component_list',2};
IgnoreList_NAME1 = [];
m=readLL(LLname,quads,IgnoreList_NAME1);
% get reference position
if ~isnumeric(reference_point)
    if ischar(reference_point)
        [stemp,optics]=readLL(LLname,reference_point,IgnoreList_NAME1);
        refpos = stemp.z_pos;
    end
else
    refpos = reference_point;
end

% sort m by z_pos
% insert code to find all magents (QUAD or SBEN) between quads{1} and reference_point
% M has to be changed to include m(:).type = 'QUAD' or 'SBEN' und
% m(:).active = 1(if in list quads) / 0 (otherwise)

twiss_goal = optics([3,2,6,5]);

n=50; N=10; K=cell(1,iter);
for ii=1:iter
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
    Kok=K2sum<=kmax^2;
    betaok=squeeze(min(min(twiss(:,[1 3],:),[],1),[],2))>=betamin&squeeze(max(max(twiss(:,[1 3],:),[],1),[],2))<=betamax;
    betxyok=squeeze(min(prod(twiss(:,[1 3],:),2),[],1))>=betaxymin;
    
    Cx=sum(squeeze(twiss(2:N+2:end,1,:).^2)'.*K{ii}.^2,2);
    Cy=sum(squeeze(twiss(2:N+2:end,3,:).^2)'.*K{ii}.^2,2);
    Cxok=Cx>=Cxmin;
    Cyok=Cy>=Cymin;
    
    allok=Kok&betaok&betxyok&Cxok&Cyok;
    K{ii}=K{ii}(allok,:);
    twiss=twiss(:,:,allok);
    twiss_store{ii} = twiss;
    
    K2sum=K2sum(allok)';
    Cx=Cx(allok)';
    Cy=Cy(allok)';
    
    [~,icxy]=sort(Cx.*Cy);
    Kmin=min(K{ii}(icxy(1:Nbest),:));
    Kmax=max(K{ii}(icxy(1:Nbest),:));
    
    best(ii,:) = icxy(1:Nbest);
    %figure(1); plot3(K{ii}(:,1),K{ii}(:,2),K{ii}(:,3),'.'); grid on; hold on; pause(0.01);
end
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