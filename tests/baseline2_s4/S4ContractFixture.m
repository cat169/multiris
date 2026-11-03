classdef S4ContractFixture < handle
    % Isolated wrapper mirror; original numerical algorithms are never edited.
    properties
        Root
        OriginalPath
        Prior
        Configuration
    end
    methods
        function obj=S4ContractFixture()
            here=fileparts(mfilename('fullpath'));repo=fileparts(fileparts(here));
            obj.OriginalPath=path;obj.Root=tempname;mkdir(obj.Root);
            inference=fullfile(repo,'baseline2/inference');
            copyfile(fullfile(inference,'infer_baseline2_frame_s4.m'),obj.Root);
            copyfile(fullfile(inference,'private'),fullfile(obj.Root,'private'));
            files=dir(fullfile(here,'doubles','*.fixture'));
            for k=1:numel(files)
                copyfile(fullfile(files(k).folder,files(k).name), ...
                    fullfile(obj.Root,erase(files(k).name,'.fixture')));
            end
            addpath(obj.Root,'-begin');
            clear infer_baseline2_frame_s4 refine_position_map estimate_channel_parameters compute_channel_fim s4_contract_state;
            obj.Prior=struct('mean',[1;1;1;2;2;2],'covariance',eye(6));
            cfg.constants.c=3e8;cfg.ofdm.subcarrier_spacing=2.5e6;
            cfg.geometry.ris_centers=zeros(3,1);
            cfg.geometry.ris_axis_x=[1;0;0];cfg.geometry.ris_axis_y=[0;1;0];
            obj.Configuration=cfg;obj.setMode('valid');
        end
        function setMode(~,mode)
            s4_contract_state('configure',mode);
        end
        function result=call(obj)
            result=infer_baseline2_frame_s4([],1,[],obj.Prior,1,obj.Configuration);
        end
        function callBadPrior(obj)
            p=obj.Prior;p.mean=p.mean.';
            infer_baseline2_frame_s4([],1,[],p,1,obj.Configuration);
        end
        function close(obj)
            path(obj.OriginalPath);
            clear infer_baseline2_frame_s4 refine_position_map estimate_channel_parameters compute_channel_fim s4_contract_state;
            assert(startsWith(obj.Root,tempdir),'s4Fixture:CleanupPath','Unexpected fixture directory.');
            if isfolder(obj.Root),rmdir(obj.Root,'s');end
        end
    end
end
