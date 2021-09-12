function phsp_out = momentum_reassign(phsp, idx, flagdist)

phsp_out = phsp;
nspawn = sum(idx);
nfix = sum(~idx);
npart = nfix+nspawn;



%Set random angular distribution
if (strcmp(flagdist,'gaussian'))   
    r0px = 5*std(phsp(~idx(:),2))
    r0py = 5*std(phsp(~idx(:),4))
    phsp_out(idx,2) = r0px*randn(nspawn,1);
    phsp_out(idx,4) = r0py*randn(nspawn,1);
elseif (strcmp(flagdist,'uniform'))
    r0px = 5*std(phsp(~idx(:),2))
    r0py = 5*std(phsp(~idx(:),4))
    r2px = rand(nspawn,1)*r0px*r0px;
    r2py = rand(nspawn,1)*r0py*r0py;
    phip = rand(nspawn,1)*2*pi;
    phsp_out(idx,2)=sqrt(r2px).*cos(phip);
    phsp_out(idx,4)=sqrt(r2py).*sin(phip);
elseif (strcmp(flagdist,'fitsim'))
    s_attempt = 4*cov(phsp(~idx(:),1:4));
    phsp_attempt = mvnrnd([0 0 0 0],s_attempt,nspawn);
%    phsp_out(idx,2) = phsp_attempt(:,2);
%    phsp_out(idx,4) = phsp_attempt(:,4);
    phsp_out(idx,1:4) = phsp_attempt(:,1:4);
elseif (strcmp(flagdist,'fitdata'))
    s_attempt = 4*cov(phsp(~idx(:),1:4));
    phsp_attempt = mvnrnd([0 0 0 0],s_attempt,nspawn);
    phsp_out(idx,[2,4]) = phsp_attempt(:,[2,4]);
elseif (strcmp(flagdist,'randomwalk'))
    r0px = 5*std(phsp(~idx(:),2))
    r0py = 5*std(phsp(~idx(:),4))
    phsp_out(idx,2) = phsp_out(idx,2)+r0px*rand(1)*cos(rand(1)*2*pi);
    phsp_out(idx,4) = phsp_out(idx,4)+r0py*rand(1)*sin(rand(1)*2*pi);
end

end
