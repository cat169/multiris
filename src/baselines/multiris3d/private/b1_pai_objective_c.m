function [cost, parts, minDistance] = b1_pai_objective_c(x, channel, geometry, c)
% Independent C objective: spatial residuals plus distance residual in metres.
% Truth fields are not passed here; parts(2) has units metres squared.
active=channel.activeRisIndices(:); residual=zeros(3*numel(active),1); minDistance=inf;
for j=1:numel(active)
    l=active(j); dT=x(1:3)-geometry.ris_centers(:,l); dR=x(4:6)-geometry.ris_centers(:,l);
    nt=norm(dT); nr=norm(dR); minDistance=min([minDistance,nt,nr]);
    if nt==0 || nr==0, cost=nan; parts=[nan nan]; return; end
    u=dT/nt+dR/nr;
    residual(3*(j-1)+(1:3))=[channel.xiEstimate(j)-geometry.ris_axis_x(:,l).'*u; channel.zetaEstimate(j)-geometry.ris_axis_y(:,l).'*u; c*channel.delayEstimate(j)-nt-nr];
end
cost=sum(residual.^2); s=reshape(residual,3,[]); parts=[sum(sum(s(1:2,:).^2)),sum(s(3,:).^2)];
end
