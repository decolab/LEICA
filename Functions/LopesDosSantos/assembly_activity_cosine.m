function [Activities] = assembly_activity_cosine(Patterns, events)

% Activities = assembly_activity(Patterns,Activitymatrix): computes the
% time course of the activity of assembly patterns defined in Patterns in
% the spike matrix Activitymatrix with single bin resolution.
% 
% Description of inputs: 	
%   Activitymatrix: spike matrix. Rows represent neurons, columns represent
%   time bins. Patterns: assembly patterns. Columns denote assembly # and
%   rows neuron #.
% 
% Description of output: 	
%   Activities: Time course of the activity of assemblies. Rows represent
%   assemblies and columns represent time bins. 
%
% This framework is described in: Lopes-dos-Santos V, Ribeiro S, Tort ABL 
% (2013) Detecting cell assemblies in large neuronal populations, Journal
% of Neuroscience Methods.
%
% Please send bug reports to vitor@neuro.ufrn.br (V�tor)

NumPatterns=size(Patterns,2);
ntime=size(events,2);
Activities = nan(NumPatterns,ntime);

% Check if any patterns detected
if NumPatterns == 0
	Activities = [];
else
	sevents=(zscore(events'))';
	for t=1:ntime
		for ass=1:NumPatterns
			Activities(ass,t)=(dot(Patterns(:,ass),sevents(:,t))/norm(Patterns(:,ass))/norm(sevents(:,t)))^2;
		end
	end
end

end