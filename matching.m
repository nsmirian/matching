function [outputArg1,outputArg2] = matching(lattice,twiss_measured,twiss_design,constraints)
%matching 
% lattice.reference_point
% lattice.quads
% lattice.Kmax
% lattice.Kmin
% lattice.k

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

if or(~exist('twiss_design','var'),isempty(twiss_design))
    twiss_goal = optics(end,[3,2,6,5]);
end

if or(~exist('constraints','var'),isempty(constraints))
    constraints = [1000,1,0,1000,0];
end

Ndrift=10;
Nmc=10;
withparfor=false;

[out(i),stability] = the_alternator(m,refpos,lattice.k,twiss_measured,twiss_goal,Ndrift,constraints,lattice.Kmax,Nmc,withparfor);

outputArg1 = inputArg1;
outputArg2 = inputArg2;
end

