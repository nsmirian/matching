function [ output,optics ] = readLL( LLname,Snames,IgnoreList_NAME1 )
%function [ output,optics ] = readLL( LLname )
% 2017__01_27 changed to strength as in long list (not length integrated!)
%             added nominal arc length (S_NOMINAL) to SBEN parameters

    [data,desc]=xlsread(LLname{1},LLname{2});
    desc_offset=size(desc)-size(data);
    header=desc(1,desc_offset(2)+1:end);
    iLENGTH  =find(strcmp(header,'LENGTH'));   %  1
    iSTRENGTH=find(strcmp(header,'STRENGTH')); %  2
    iS       =find(strcmp(header,'S'));        %  6
    iX       =find(strcmp(header,'X'));        %  8
    iY       =find(strcmp(header,'Y'));        %  9
    iZ       =find(strcmp(header,'Z'));        % 10
    iALFX    =find(strcmp(header,'ALFX'));
    iALFY    =find(strcmp(header,'ALFY'));
    iBETX    =find(strcmp(header,'BETX'));     % 21
    iBETY    =find(strcmp(header,'BETY'));     % 24
    iMUX     =find(strcmp(header,'MUX'));
    iMUY     =find(strcmp(header,'MUY'));
    % get required sections -----------------------------------------------
    if iscell(Snames)
        ilines=find(strcmp(desc(:,3),Snames{1})|strcmp(desc(:,4),Snames{1}));
        for i=2:size(Snames,2)
            ilines=[ilines;find(strcmp(desc(:,3),Snames{i})|strcmp(desc(:,4),Snames{i}))];
        end
    else
        ilines=find(strcmp(desc(:,3),Snames)|strcmp(desc(:,4),Snames));
    end
    
    [~,idx]=sort(data(ilines-desc_offset(1),iS));
    ilines=ilines(idx);
    data=data(ilines-desc_offset(1),:);
    NAME1=desc(ilines,strcmp(desc(1,:),'NAME1'));
    clear desc header
    % optics --------------------------------------------------------------
    optics=data(:,[iS iALFX iBETX iMUX iALFY iBETY iMUY]);
    % list of ignored elements
    ilines=false(size(NAME1));
    for i=1:size(IgnoreList_NAME1,2)
        ilines=ilines | strcmp(NAME1,IgnoreList_NAME1{i});
    end
    Nlines=size(ilines,1);
    %output=cell(Nlines,1);
    for i=1:Nlines
        output(i)=struct('effective_length',data(i,iLENGTH),...
                         'strength_sp',data(i,iSTRENGTH)/data(i,iLENGTH),...
                         'z_pos',data(i,iZ));

    end
end
