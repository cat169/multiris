function steering = upa_steering(xi, zeta, cfg)
%{
根据 xi、zeta 和阵元间距生成零起点、未归一化的 UPA 导向矢量。输出使用
kron(x 轴响应,y 轴响应)。
%}

indexX = (0:(cfg.ris.nx - 1)).';
indexY = (0:(cfg.ris.ny - 1)).';

steeringX = exp(1i * (2 * pi / cfg.constants.lambda) * ...
    cfg.ris.dx * indexX * xi);
steeringY = exp(1i * (2 * pi / cfg.constants.lambda) * ...
    cfg.ris.dy * indexY * zeta);

steering = kron(steeringX, steeringY);
end
