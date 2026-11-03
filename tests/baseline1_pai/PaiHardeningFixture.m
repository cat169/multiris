classdef PaiHardeningFixture < handle
    % Isolated source mirror and solver doubles; production files are read only.
    properties
        Root
        OriginalPath
        Configuration
        Channel
    end
    methods
        function obj=PaiHardeningFixture()
            here=fileparts(mfilename('fullpath')); repo=fileparts(fileparts(here));
            obj.OriginalPath=path; obj.Root=tempname; mkdir(obj.Root);
            source=fullfile(repo,'src','baselines','multiris3d');
            copyfile(fullfile(source,'*.m'),obj.Root); mkdir(fullfile(obj.Root,'private'));
            copyfile(fullfile(source,'private','*.m'),fullfile(obj.Root,'private'));
            fixtures=dir(fullfile(here,'doubles','*.fixture'));
            for k=1:numel(fixtures)
                name=erase(fixtures(k).name,'.fixture'); copyfile(fullfile(fixtures(k).folder,fixtures(k).name),fullfile(obj.Root,name));
            end
            addpath(obj.Root,'-begin'); clear infer_baseline1_frame_projection_assisted estimate_position_coarse refine_position_mle_eq57_identity estimate_channel_parameters pai_hardening_stub_state;
            cfg.constants.c=3e8; cfg.ris.count=2;
            cfg.geometry.ris_centers=[0 0;0 2;0 0]; cfg.geometry.ris_axis_x=[1 1;0 0;0 0]; cfg.geometry.ris_axis_y=[0 0;0 0;1 1]; cfg.geometry.ris_axis_normal=[0 0;1 1;0 0];
            cfg.geometry.static_user=zeros(3,1); cfg.geometry.mobile_user=zeros(3,1);
            cfg.paperBaselineA=struct('oracleSupport',[1;2],'lowerBound',[1;0;0],'upperBound',[2;0;0],'gridStep',1);
            obj.Configuration=cfg; obj.Channel=struct('activeRisIndices',[1;2],'delayEstimate',[3e-8;4e-8],'xiEstimate',[1.5;-.4],'zetaEstimate',[.5;.1],'separatedNoiseVariance',1);
            pai_hardening_stub_state('configure',struct('mode','finite','channel',obj.Channel));
        end
        function setMode(obj,mode)
            pai_hardening_stub_state('configure',struct('mode',mode,'channel',obj.Channel));
        end
        function out=call(obj)
            out=infer_baseline1_frame_projection_assisted([],1,[],obj.Configuration,obj.Channel);
        end
        function out=callDefault(obj)
            out=infer_baseline1_frame_projection_assisted([],1,[],obj.Configuration);
        end
        function value=state(~)
            value=pai_hardening_stub_state('get',[]);
        end
        function close(obj)
            path(obj.OriginalPath); clear infer_baseline1_frame_projection_assisted estimate_position_coarse refine_position_mle_eq57_identity estimate_channel_parameters pai_hardening_stub_state;
            if isfolder(obj.Root), rmdir(obj.Root,'s'); end
        end
    end
end
