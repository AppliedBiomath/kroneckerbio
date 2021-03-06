function con = Uzero(m)
%UZERO Default structure for Kronecker Bio experimental conditions
% 
%   This function returns a valid, empty, experiment. It is used for
%   making experiments constructed by different means consistent so that
%   their structures can be combined.
%
%   Inputs
%   m: [ Kronecker model structure scalar ]
%
%   Outputs:
%   con: [ Kronecker experimental condition structure scalar ]
%       The experiment returned by this function is pretty useless, but the
%       fields can be overwritten. Their meanings are given below.
%       .Type 'Experiment'
%       .Name [ string ]
%           An aribitrary name for the experiment
%       .tF [ nonnegative scalar ]
%           The time at which this experiment ends.
%       .s [ nonnegative vector ns ]
%           The values of the seed parameters
%       .u [ handle @(t) returns nonnegative matrix nu by numel(t) ]
%           A function handle that returns the value of the inputs at any
%           time t. This function should accept t as a vector.
%       .q [ real vector nq ]
%           The values of the input control parameters
%       .dudq [ handle @(t) returns real matrix nu by nq]
%           A function handle that returns the derivative of each input
%           with respect to each input control parameter
%       .nq [ whole scalar ]
%           The number of input control parameters
%       .nqu [ whole vector nv ]
%           The number of input control parameters associated with each
%           input
%       .SteadyState [ logical scalar ]
%           Declares if the system should be run to steady state before
%           the experiment begins
%       .Periodic [ logical scalar ]
%           Declares if the system should be run to a periodic steady state
%           for the duration of the experiment
%       .Discontinuities [ nonegative vector ]
%           Any discontinuous times in the u function must be listed here
%           in order to ensure sucessful integration of the system.
%       .Update [handle @(s,q) returns struct scalar ]
%           This function handle allows the parameter values of the
%           experiment to be changed. Each vector x0 and q must have the
%           same size as their appropriate counterparts. The structure
%           returned is identical to the one to which Update is attached
%           except that the parameter values have been changed and the
%           appropriate matrices and function handles have been updated.

% (c) 2011 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

% Special case of numeric inputs
if isnumeric(m)
    con = emptystruct(m, 'Type', 'Name', 'tF', 's', 'u', 'q', 'dudq', 'nq', 'SteadyState', 'Periodic', 'Discontinuities', 'Update');
    return
end

% Constants
nu = m.nu;
s = m.s;

% Gut m
temp = m;
clear m
m.nu = temp.nu;
m.ns = temp.ns;
m.s = temp.s;
clear temp

con.Type = 'Experiment';
con.Name = 'UnnamedExperiment';
con.tF = 0;
con.s = s;
con.u  = @(t)zeros(nu,1);
con.q  = zeros(0,1);
con.dudq = @(t)zeros(nu,0);
con.nq = 0;
con.SteadyState = false;
con.Periodic = false;
con.Discontinuities = zeros(0,1);
con.Update = @Update;

    function con = Update(s, q)
        assert(numel(q) == 0, 'KroneckerBio:Experiment:Update:InvalidInputControl', 'An input control parameter vector was supplied with %i elements, but the input of this experiment is controlled by 0 input parameters. They must be the same.', numel(q))
        con = Experiment(m, 0, s, false, false, zeros(nu,1), [], q);
    end
end