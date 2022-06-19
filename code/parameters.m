function params = parameters(varargin)
    if strcmp(varargin{1}, 'p')
        params = varargin{2};
    else
        params.dt = 5e-5;
        params.H1 = 1; params.L1 = -1; 
        params.H2 = 1; params.L2 = -1; % set thresholds
        params.mu1 = 1; params.mu2 = -1; % set drift rates
        params.D = 1; % set diffusion coefficient
        params.qp1 = 1; % set positive-threshold "kick" strength
        params.qp2 = 1;
        params.qn1 = 1; % set negative-threshold "kick" strength
        params.qn2 = 1;
        params.R1p = 1; params.R1n = 1; % set agent 1 rewards at each threshold
        params.R2p = 1; params.R2n = 1; % set agent 2 rewards at each threshold
    end

    if any(strcmp(varargin, 'dt'))
        params.dt = varargin{find(strcmp(varargin,'dt'))+1};
    end
    if any(strcmp(varargin, 'H1'))
        params.H1 = varargin{find(strcmp(varargin,'H1'))+1};
    end
    if any(strcmp(varargin, 'L1'))
        params.L1 = varargin{find(strcmp(varargin,'L1'))+1};
    end
    if any(strcmp(varargin, 'H2'))
        params.H2 = varargin{find(strcmp(varargin,'H2'))+1};
    end
    if any(strcmp(varargin, 'L2'))
        params.L2 = varargin{find(strcmp(varargin,'L2'))+1};
    end
    if any(strcmp(varargin, 'mu1'))
        params.mu1 = varargin{find(strcmp(varargin,'mu1'))+1};
    end
    if any(strcmp(varargin, 'mu2'))
        params.mu2 = varargin{find(strcmp(varargin,'mu2'))+1};
    end
    if any(strcmp(varargin, 'qp1'))
        params.qp1 = varargin{find(strcmp(varargin,'qp1'))+1};
    end
    if any(strcmp(varargin, 'qp2'))
        params.qp2 = varargin{find(strcmp(varargin,'qp2'))+1};
    end
    if any(strcmp(varargin, 'qn1'))
        params.qn1 = varargin{find(strcmp(varargin,'qn1'))+1};
    end
    if any(strcmp(varargin, 'qn2'))
        params.qn2 = varargin{find(strcmp(varargin,'qn2'))+1};
    end
    if any(strcmp(varargin, 'R1p'))
        params.R1p = varargin{find(strcmp(varargin,'R1p'))+1};
    end
    if any(strcmp(varargin, 'R1n'))
        params.R1n = varargin{find(strcmp(varargin,'R1n'))+1};
    end
    if any(strcmp(varargin, 'R2p'))
        params.R2p = varargin{find(strcmp(varargin,'R2p'))+1};
    end
    if any(strcmp(varargin, 'R2n'))
        params.R2n = varargin{find(strcmp(varargin,'R2n'))+1};
    end
end