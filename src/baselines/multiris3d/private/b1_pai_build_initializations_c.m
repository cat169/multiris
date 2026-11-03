function [initializations,diag]=b1_pai_build_initializations_c(channel,geometry,search,c)
grid=b1_pai_grid(search); active=channel.activeRisIndices(:); pairs=nchoosek(1:numel(active),2); n=size(grid,2); L=numel(active); P=size(pairs,1);
diag=struct('directionNormBefore',nan(L,n),'directionNormAfter',nan(L,n),'projectionDistance',nan(L,n),'normalComponent',nan(L,n),'projected',false(L,n),'pairRank',nan(P,n),'pairCondition',nan(P,n),'pairClosestDistance',nan(P,n),'pairRanges',nan(2*P,n),'nonfiniteIntersection',false(1,n),'invalidObjective',false(1,n),'txAtRisCenter',false(1,n),'gridCount',n);
template=struct('x0',nan(6,1),'gridIndex',nan,'risIndex',0,'remainingDistance',nan,'initialObjective',nan,'projectionError',nan,'maxPairCondition',nan,'delayL1Metres',nan);
pool=repmat(template,1,n); count=0;
for gi=1:n
    tx=grid(:,gi); dirs=zeros(3,L);
    for j=1:L
        l=active(j); d=tx-geometry.ris_centers(:,l); distance=norm(d);
        if distance==0, diag.txAtRisCenter(gi)=true; break; end
        u=d/distance;
        ab=[channel.xiEstimate(j)-geometry.ris_axis_x(:,l).'*u; channel.zetaEstimate(j)-geometry.ris_axis_y(:,l).'*u];
        old=ab; radius=norm(ab); if radius>1, ab=ab/radius; end
        normal=sqrt(max(0,1-sum(ab.^2)));
        dirs(:,j)=ab(1)*geometry.ris_axis_x(:,l)+ab(2)*geometry.ris_axis_y(:,l)+normal*geometry.ris_axis_normal(:,l);
        diag.directionNormBefore(j,gi)=radius; diag.directionNormAfter(j,gi)=norm(ab); diag.projectionDistance(j,gi)=norm(ab-old); diag.normalComponent(j,gi)=normal; diag.projected(j,gi)=radius>1;
    end
    if diag.txAtRisCenter(gi), continue; end
    assert(all(diag.directionNormAfter(:,gi)<=1+1e-12),'fallback:ProjectionConstraint','Projection outside unit disk.');
    points=nan(3,P); valid=true;
    for p=1:P
        i=pairs(p,1); j=pairs(p,2); A=[dirs(:,i),-dirs(:,j)]; diag.pairRank(p,gi)=rank(A); diag.pairCondition(p,gi)=cond(A'*A);
        if diag.pairRank(p,gi)<2, valid=false; continue; end
        ri=geometry.ris_centers(:,active(i)); rj=geometry.ris_centers(:,active(j)); ranges=A\(rj-ri);
        pointI=ri+ranges(1)*dirs(:,i); pointJ=rj+ranges(2)*dirs(:,j); points(:,p)=pointI;
        diag.pairRanges(2*(p-1)+(1:2),gi)=ranges; diag.pairClosestDistance(p,gi)=norm(pointI-pointJ);
        if any(~isfinite([ranges;pointI;pointJ])), valid=false; diag.nonfiniteIntersection(gi)=true; end
    end
    if ~valid, continue; end
    rx=mean(points,2); x0=[tx;rx]; [cost,~,~]=b1_pai_objective_c(x0,channel,geometry,c);
    if ~isreal(x0) || any(~isfinite(x0)) || ~isfinite(cost), diag.invalidObjective(gi)=true; continue; end
    count=count+1; pool(count)=template; pool(count).x0=x0; pool(count).gridIndex=gi; pool(count).initialObjective=cost;
    pool(count).projectionError=sum(diag.projectionDistance(:,gi).^2); pool(count).maxPairCondition=max(diag.pairCondition(:,gi));
    centers=geometry.ris_centers(:,active); delayResidual=vecnorm(tx-centers).'+vecnorm(rx-centers).'-c*channel.delayEstimate;
    pool(count).delayL1Metres=sum(abs(delayResidual));
end
diag.validInitializations=count; initializations=b1_pai_rank_initializations(pool(1:count));
end
