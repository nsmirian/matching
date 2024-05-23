function out = the_alternator(m,refpos,k_rand,twiss_in,twiss_goal,N,zlist)  %,twiss_in,twiss_goal
%  out = the_alternator(quads,reference_point,twiss_in,twiss_goal)
%  returns values for k's of quads to achive matching
%   input:
%                quads     - cell vetor of quadrupole names (it is assumed to have only drifts between them)
%                reference point  - start point for matching either a string name or a position in [m] (e.g. fodo measurement reference point)
%                twiss_in    -  measured twiss parameters (beta_x, alpha_x, beta_y, alpha_y)
%                twiss_goal  -  matching target parameters

    %for ii=1:length(m)
    %    k_orig(ii) = m(ii).strength_sp;
    %end

    k_orig = [ m( [m(:).active]==1 ).strength_sp ];
    
    matrix_orig = assemble_matrix(m,k_orig,refpos);
    %out = get_matching_matrix(k_orig*1.001,q_leff,D)
    twiss_back = track(matrix_orig{end}^(-1),twiss_in);
    twiss_orig = track(matrix_orig,twiss_back);
    
    options=optimset('MaxFunEvals',length(m)*1000,'MaxIter',length(m)*1000,'TolFun',5e-4);
    k_new = fminsearch(@(k) match_cond(twiss_back,twiss_goal,get_matching_matrix_curly(m,k,refpos,zlist),zlist),k_rand,options);

    mismatch_xy = match_cond(twiss_back,twiss_goal,get_matching_matrix_curly(m,k_new,refpos));

    [outtmp,z,quadindx] = assemble_matrix(m,k_new,refpos);
    matched=track(outtmp , twiss_back) ;

%     figure(24)
%     %subplot(2,1,1)
%     plot(z,twiss_orig(:,1),'b')
%     hold on
%     plot(z,twiss_orig(:,3),'r')
%     plot(z,matched(:,1),'b--')
%     plot(z,matched(:,3),'r--')
%     %hold off
%     %subplot(2,1,2)
%     %bar([k_orig;k_new]');
%     refresh 

    out.k_new = k_new;
    out.matched = matched;
    out.twiss_orig = twiss_orig;
    out.k_orig = k_orig;
    out.mismatch_xy = mismatch_xy;
    out.z = z;
    out.quadindx = quadindx;

    function [out,quadindx] = assemble_matrix_simple(m,k,refpos)
        % generate matching section
        k_temp = [m(:).strength_sp];
        indx = [m(:).active]==1;
        k_temp(indx) = k;
        
        alpha = 0;
        
        quadindx = false;
        out{1} = eye(6);
        count = 2;
        for i=1:length(m)
            switch m(i).type
                case 'QUAD'
                    M = Quadrupole(k_temp(i),m(i).effective_length);
                    out{count} = M*out{count-1};count=count+1;
                    quadindx = [quadindx true];
                case 'SBEN'
                    Dalpha = m(i).strength_sp*m(i).effective_length;
                    alpha = alpha +Dalpha;
                    r = m(i).effective_length/sin(Dalpha);
                    temp = Poleface(Dalpha*180/pi,r) * Sector(Dalpha*180/pi,r) * Poleface(Dalpha*180/pi,r) ;
                    % exchange x and y to describe vertical bends
                    M = temp;
                    M(1:2,1:2) = temp(3:4,3:4);
                    M(3:4,3:4) = temp(1:2,1:2);
                    out{count} = M*out{count-1};count=count+1;
                    quadindx = [quadindx false];
            end
            if i<length(m)
                D = Drift( (m(i+1).z_pos - m(i).z_pos -0.5*(m(i).effective_length + m(i+1).effective_length)) /cos(alpha) );
            else
                D = Drift( (refpos - m(end).z_pos + m(end).effective_length/2)/cos(alpha) );
            end
            out{count} = D*out{count-1};count =count+1;
            quadindx = [quadindx false];
        end
    end
    
    function [out,z,quadindx] = assemble_matrix(m,k,refpos)
        % generate matching section
        k_temp = [m(:).strength_sp]; 
        indx = [m(:).active]==1;
        k_temp(indx) = k;
        
        %out=cell(1,length(m));
        z=m(1).z_pos-m(1).effective_length/2;
        
        alpha = 0;
        
        quadindx = false;
        out{1} = eye(6);
        count = 2;
        for i=1:length(m)           
            switch m(i).type
                case 'QUAD'
                    M = Quadrupole(k_temp(i),m(i).effective_length/2);
                    out{count} = M*out{count-1};count=count+1;
                    z = [z, z(end)+m(i).effective_length/2];
                    out{count} = M*out{count-1};count=count+1;
                    z = [z, z(end)+m(i).effective_length/2];
                    quadindx = [quadindx true false];          
                case 'SBEN'
                    Dalpha = m(i).strength_sp*m(i).effective_length;
                    alpha = alpha +Dalpha;
                    r = m(i).effective_length/sin(Dalpha);
                    temp = Poleface(Dalpha*180/pi,r) * Sector(Dalpha*180/pi,r) * Poleface(Dalpha*180/pi,r) ;
                    % exchange x and y to describe vertical bends
                    M = temp;
                    M(1:2,1:2) = temp(3:4,3:4);
                    M(3:4,3:4) = temp(1:2,1:2);
                    out{count} = M*out{count-1};count=count+1;
                    quadindx = [quadindx false];
                    z = [z, z(end)+m(i).effective_length];
            end
            if i<length(m)
               D = Drift( (m(i+1).z_pos - m(i).z_pos -0.5*(m(i).effective_length + m(i+1).effective_length))/N /cos(alpha) );
            else
               D = Drift( (refpos - m(end).z_pos + m(end).effective_length/2)/N/cos(alpha) );
            end
            for k=1:N
                out{count} = D*out{count-1};count =count+1;
                quadindx = [quadindx false];
                z = [z, z(end)+D(1,2)*cos(alpha)];
            end
        end
    end

%     function Q = make_quads(q_leff,k)
%         for i=1:length(q_leff)
%             Q{i} = Quadrupole(k(i),q_leff(i));
%         end
%     end

    function out = get_matching_matrix(m,k,refpos,zlist)
        [outtemp,ztemp] = assemble_matrix(m,k,refpos);
        
        if nargin==4
            zcomp = @(x,y)abs(x-y);
            % finde outtemp um die punkte zlist
            [~,indx] = min(bsxfun(zcomp,zlist(:,1),ztemp),[],2);
            %indx = [indx indx+1];
            %indx = reshape(indx',1,2*size(zlist,1));
            out = outtemp([indx(:);length(outtemp)]);
        else
            out = outtemp{end};
        end
    end

    function out = get_matching_matrix_curly(k,leff,D,zlist)
        % aendere funtion so dass cell von matrizen gecurled wird
        if nargin==4
            temp = get_matching_matrix(k,leff,D,zlist);
            out = [];
            for i=1:length(temp)
                out = cat(3,out,CurlyM(temp{i}));
            end
        else
            temp = get_matching_matrix(k,leff,D);
            out = CurlyM(temp);
        end
    end
end

function out = track(M,twissin)
% twiss in (beta_x,alpha_x,beta_y,alpha_y)
    twissin=twissin(:)';
    twiss    = [twissin(1:2)  , (1+twissin(2)^2)/twissin(1)     ,twissin(3:4)  , (1+twissin(4)^2)/twissin(3) ];
    if iscell(M)
        for i=1:length(M)
            twissout = CurlyM(M{i})*twiss(:);
            out(i,:) = twissout([1:2,4:5]);
        end
    else
        twissout = CurlyM(M)*twiss(:);
        out = twissout([1:2,4:5])';
    end
end

function penalty = match_cond(twissin,twiss_out,mat,zlist)
% twiss in (beta_x,alpha_x,beta_y,alpha_y)
    twissin=twissin(:)';
    twiss    = [twissin(1:2)  , (1+twissin(2)^2)/twissin(1)     ,twissin(3:4)  , (1+twissin(4)^2)/twissin(3) ];
    twissout = [twiss_out(1:2), (1+twiss_out(2)^2)/twiss_out(1) ,twiss_out(3:4), (1+twiss_out(4)^2)/twiss_out(3) ];

    % aendere so dass temp fuer alle matrizen asgerechnet wird
    if nargin==4
        nm=size(mat,3);
        mat=reshape(permute(mat,[3 1 2]),6*nm,6);
        temp=reshape(mat*twiss(:),nm,6);
        mismatch_x = 0.5 * (twissout(1)*temp(end,3)   -  2*twissout(2)*temp(end,2) + twissout(3)*temp(end,1) );
        mismatch_y = 0.5 * (twissout(4)*temp(end,6)   -  2*twissout(5)*temp(end,5) + twissout(6)*temp(end,4) );
        penalty = mismatch_x*mismatch_y  + optim_function(diag(temp(1:end-1,zlist(:,2))),zlist(:,3),zlist(:,4))'*zlist(:,5);
    else
        temp = mat*twiss(:);
        mismatch_x = 0.5 * (twissout(1)*temp(3)   -  2*twissout(2)*temp(2) + twissout(3)*temp(1) );
        mismatch_y = 0.5 * (twissout(4)*temp(6)   -  2*twissout(5)*temp(5) + twissout(6)*temp(4) );
        penalty = mismatch_x*mismatch_y;
    end

    % symbolisch (i index in zlist zeile) + \sum_i alpha_i*optim_function( temp(i,zlist(i,2)) ,zlist(i,3),zlist(i,4)) 
end

function f = optim_function(x,minx,maxx)
    f=(x<minx).*(sqrt((x-minx).^2+9)-3)+(x>maxx).*(sqrt((x-maxx).^2+9)-3);
end
