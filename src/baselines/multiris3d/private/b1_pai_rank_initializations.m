function selected=b1_pai_rank_initializations(pool)
% Numeric sorting and position deduplication only. No truth input.
selected=pool([]); if isempty(pool), return; end
key=[[pool.initialObjective].',[pool.gridIndex].',[pool.risIndex].']; [~,order]=sortrows(key,[1 2 3]);
for k=order.'
    candidate=pool(k);
    if ~isempty(selected) && any(vecnorm([selected.x0]-candidate.x0)<=1e-6), continue; end
    selected(end+1)=candidate; %#ok<AGROW>
    if numel(selected)==8, break; end
end
end
