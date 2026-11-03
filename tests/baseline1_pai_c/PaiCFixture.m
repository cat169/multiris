classdef PaiCFixture
    % Synthetic channels only; no MLE, observation generation or MC.
    methods (Static)
        function [cfg,ch,x]=data(pattern)
            cfg.constants.c=3e8;cfg.ris.count=3;
            cfg.geometry.ris_centers=[-5 5 0;0 0 5;0 0 0];
            cfg.geometry.ris_axis_x=repmat([1;0;0],1,3);
            cfg.geometry.ris_axis_y=repmat([0;1;0],1,3);
            cfg.geometry.ris_axis_normal=repmat([0;0;1],1,3);
            cfg.geometry.static_user=[0;0;4];cfg.geometry.mobile_user=[1;2;3];
            active=find(char(pattern).'-'0');
            cfg.paperBaselineA=struct('oracleSupport',active,'lowerBound',[-.2;-.2;3.8], ...
                'upperBound',[.2;.2;4.2],'gridStep',.2);
            x=[cfg.geometry.static_user;cfg.geometry.mobile_user];
            ch=struct('activeRisIndices',active,'delayEstimate',zeros(numel(active),1), ...
                'xiEstimate',zeros(numel(active),1),'zetaEstimate',zeros(numel(active),1), ...
                'separatedNoiseVariance',1,'exitflag',ones(numel(active),1));
            residual=PaiCFixture.secondsResidual(x,ch,cfg);
            ch.xiEstimate=-residual(1:3:end)+.001;
            ch.zetaEstimate=-residual(2:3:end)-.001;
            ch.delayEstimate=-residual(3:3:end)+1e-11;
        end
        function residual=secondsResidual(x,ch,cfg)
            residual=zeros(3*numel(ch.activeRisIndices),1);
            for j=1:numel(ch.activeRisIndices)
                l=ch.activeRisIndices(j);a=x(1:3)-cfg.geometry.ris_centers(:,l);b=x(4:6)-cfg.geometry.ris_centers(:,l);
                direction=a/norm(a)+b/norm(b);
                residual(3*(j-1)+(1:3))=[ch.xiEstimate(j)-cfg.geometry.ris_axis_x(:,l).'*direction; ...
                    ch.zetaEstimate(j)-cfg.geometry.ris_axis_y(:,l).'*direction; ...
                    ch.delayEstimate(j)-(norm(a)+norm(b))/cfg.constants.c];
            end
        end
        function cost=independentCost(x,ch,cfg)
            residual=PaiCFixture.secondsResidual(x,ch,cfg);residual(3:3:end)=cfg.constants.c*residual(3:3:end);
            cost=sum(residual.^2);
        end
        function [absolute,relative]=finiteDifference(pattern,near)
            [cfg,ch,x]=PaiCFixture.data(pattern);
            if near,x(1:3)=cfg.geometry.ris_centers(:,1)+[1;2;3]*1e-10;end
            [~,J]=baseline1_distance_scaled_residual(x,ch.delayEstimate,ch.xiEstimate,ch.zetaEstimate,ch.activeRisIndices,cfg);
            numeric=zeros(size(J));
            for j=1:6
                step=1e-5*max(1,abs(x(j)));if near,step=1e-14;end
                plus=x;minus=x;plus(j)=plus(j)+step;minus(j)=minus(j)-step;
                rp=baseline1_distance_scaled_residual(plus,ch.delayEstimate,ch.xiEstimate,ch.zetaEstimate,ch.activeRisIndices,cfg);
                rm=baseline1_distance_scaled_residual(minus,ch.delayEstimate,ch.xiEstimate,ch.zetaEstimate,ch.activeRisIndices,cfg);
                numeric(:,j)=(rp-rm)/(2*step);
            end
            absolute=max(abs(J(:)-numeric(:)));relative=norm(J-numeric,'fro')/max(1,norm(J,'fro'));
            fprintf('PAIC_JACOBIAN pattern=%s near=%d maxAbs=%.12g relative=%.12g\n',pattern,near,absolute,relative);
        end
        function result=entry(fallback)
            [cfg,ch,~]=PaiCFixture.data("110");
            if fallback,ch.xiEstimate=[3;2.8];ch.zetaEstimate=[.5;-.7];end
            result=infer_baseline1_frame_pai_c([],1,[],cfg,ch);
        end
        function [initialDifference,finalDifference,sorted,selected]=candidateChecks(result)
            [cfg,~,~]=PaiCFixture.data("110");ch=result.channel;candidates=result.candidates;
            initialDifference=zeros(1,numel(candidates));finalDifference=initialDifference;
            for j=1:numel(candidates)
                q=candidates(j);expected=PaiCFixture.independentCost(q.initialization.x0,ch,cfg);
                initialDifference(j)=abs(q.initialObjective-expected)/max(1,abs(expected));
                if q.valid
                    cost=PaiCFixture.independentCost(q.estimate,ch,cfg);
                    finalDifference(j)=max(abs([q.finalObjective-cost,q.formal.finalCost-cost,q.formal.initialCost-expected]))/max(1,abs(expected));
                end
            end
            sorted=issorted([candidates.initialObjective]);valid=find([candidates.valid]);
            selected=false;if ~isempty(valid),[~,j]=min([candidates(valid).finalObjective]);selected=result.selectedIndex==valid(j);end
        end
        function valid=protectedHashes()
            repo=fileparts(fileparts(fileparts(mfilename('fullpath'))));
            entries=jsondecode(fileread(fullfile(repo,'tests/baseline1_pai_c/protected_sources.json')));
            valid=true;
            for j=1:numel(entries)
                bytes=unicode2native(strrep(fileread(fullfile(repo,entries(j).path)),sprintf('\r\n'),sprintf('\n')),'UTF-8');
                digest=java.security.MessageDigest.getInstance('SHA-256');digest.update(bytes);
                hex=lower(reshape(dec2hex(typecast(digest.digest(),'uint8'),2).',1,[]));
                valid=valid&&strcmp(hex,entries(j).sha256);
            end
        end
    end
end
