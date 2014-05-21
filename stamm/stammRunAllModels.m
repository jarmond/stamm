function stammRunAllModels(models,mode,varargin)
% RUNALLMODELS Runs fit automatically on all models
%
%   RUNALLMODELS(MODELS,MODE,...) Runs fit automatically on all MODELS. MODE
%   specifies fitting model and can be one of 'perturb', 'timepointko',
%   'normal', 'multicluster' (see README.txt). MODELS is a cell list of
%   strings. The appropriate parameters to the standard function should be
%   supplied following MODE, excluding the model parameter.

switch mode
  case 'perturb'
    fn=@call_perturbation;
  case 'timepointko'
    fn=@call_timepointKnockout;
  case 'normal'
    fn=@call_fitCluster;
  case 'multicluster'
    fn=@call_multicluster;
  otherwise
    error('Unknown mode')
end


for i=1:length(models)
    fprintf('Running model %s\n',models{i});
    fn(models{i});
end


    function call_perturbation(model)
    stammPerturbation(varargin{1},model,varargin{2},varargin{3}, ...
                     varargin{4},varargin{5},varargin{6},varargin{7});
        
    end

    function call_timepointKnockout(model)
    stammTimepointKnockout(varargin{1},model,varargin{2},varargin{3}, ...
                          varargin{4},varargin{5});
    end

    function call_fitCluster(model)
    stammFitCluster(varargin{1},model,varargin{2},varargin{3},varargin{4},varargin{5});
    end

    function call_multicluster(model)
        for i=1:varargin{3}
            fprintf('Fitting clustering %d of %d\n',i,varargin{3});
            cfile=num2str(i,'/multicluster_%02d.mat');
            stammClusterGenes(varargin{1},varargin{2},[varargin{5} cfile]);
            stammFitCluster(varargin{1},model,varargin{4},[varargin{5} cfile], ...
                                    varargin{5},varargin{6});
        end
    end   
end
